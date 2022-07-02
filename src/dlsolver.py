#!/usr/bin/env python
##
## dlsolver.py
##
##  Created on: Jun 18, 2020
##      Author: Jinqiang Yu
##      E-mail: jinqiang.yu@monash.edu
##

#
#==============================================================================
from pysat.solvers import Glucose3
from pysat.formula import CNF, WCNF
from pysat.examples.rc2 import RC2, RC2Stratified
import resource
from pysat.card import *
import collections
import decimal

#
#==============================================================================
class DLSolver():
    """
        Class implementing the MaxSAT approach to generating decision lists
    """
    def __init__(self, data, options):
        """
            Constructor.
        """
        self.init_stime = resource.getrusage(resource.RUSAGE_SELF).ru_utime
        self.init_ctime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime

        self.data = data
        self.options = options

        # init variable id pool
        self.reset_idpool()

        # samples divided into classes
        self.samps = {self.data.fvmap.dir[(self.data.names[-1], v)]: [] for v in sorted(self.data.feats[-1])}

        # covers by class
        self.covrs = []

        for i, s in enumerate(self.data.samps):
            self.samps[s[-1]].append(i)

        # binarize dataset if necessary
        self.binarize()
        # get missing values for each sample
        self.get_missing()

        self.cost = 0

    def reset_idpool(self):
        """
            Reset the pool of variable ids.
        """

        self.idpool = IDPool(start_from=1)

    def binarize(self):
        """
            Do one-hot encoding.
        """

        FFMap = collections.namedtuple('FFMap', ['dir', 'opp'])
        self.ffmap = FFMap(dir={}, opp={})

        curr_id = 0
        vfmap = {}  # mapping from a feature id to a list of feature ids
        for r, (name, feats) in enumerate(zip(self.data.names[:-1], self.data.feats[:-1])):
            fgroup = []
            if len(feats) != 2:
                vars_ = sorted([self.data.fvmap.dir[name, v] for v in feats])
                for i, var in enumerate(vars_):
                    vfmap[var] = [-v for v in vars_]
                    vfmap[var][i] = var

                    self.ffmap.opp[i + curr_id] = var
                    fgroup.append(i + curr_id)

                curr_id += len(feats)
            else:
                var = self.data.fvmap.dir[name, sorted(feats)[0]]
                vfmap[var] = [var]
                vfmap[-var] = [-var]
                self.ffmap.opp[curr_id] = var
                fgroup.append(curr_id)

                curr_id += 1

            self.ffmap.dir[r] = fgroup

        # rewriting samples
        for i in range(len(self.data.samps)):
            samp, out = self.data.samps[i][:-1], self.data.samps[i][-1]

            self.data.samps[i] = []
            for l in samp:
                self.data.samps[i].extend(vfmap[l])

            self.data.samps[i].append(out)

        self.nof_feats = curr_id



    def get_missing(self):
        """
            Get a list of missing values for each sample.
        """

        self.data.vmiss = []

        for s in self.data.samps:
            missing = []

            if len(s) < self.nof_feats + 1:
                r = i = 0
                while i < len(s) - 1:
                    if r in self.ffmap.dir[self.data.nm2id[self.data.fvmap.opp[abs(s[i])][0]]]:
                        i += 1
                    else:
                        missing.append(r)

                    r += 1

                # adding the rest of the features
                missing.extend(range(r, self.nof_feats))

            # set is needed for testing inclusion
            self.data.vmiss.append(set(missing))

    def ordercls(self):
        if self.options.mode == 'sep':
            #if self.options.clsorder in ('acry', 'cost'):
            if self.options.clsorder != 'maj':
                self.options.clsorder = 'maj'
                print('only support sorting classes by major classes')

            if self.options.clsorder == 'maj':
                self.labels = sorted(self.labels,
                                     key=lambda l: len(self.samps[l]),
                                     reverse=self.options.clsdown)

        else:
            self.labels = [0]

    def compute(self):
        #print('self.nof_feats:', self.nof_feats)
        #print('self.data.samps:', self.data.samps)
        #print('self.samps:', self.samps)
        #print('\nself.data.fvmap:',self.data.fvmap)
        #print('\nself.ffmap:',self.ffmap)
        self.labels = self.samps.keys()
        self.nof_labls = len(self.samps.keys())
        self.nof_labls = self.nof_labls if self.nof_labls > 2 else 1
        self.samps_ = {label: set(self.samps[label]) for label in self.samps}
        solver = RC2Stratified if self.options.approach == 'sparse' else RC2
        self.ordercls()

        self.time = 0.0
        self.nof_misses = 0

        for label in self.labels:
            self.label = label
            if self.nof_labls == 1:
                self.potential_labels = set([self.nof_feats+1])
            elif self.options.mode != 'sep':
                self.potential_labels = set(range(self.nof_feats + 1, self.nof_feats + self.nof_labls + 1))
            else:
                self.potential_labels = set([label])
            self.potential_feats = sorted(list(range(1, self.nof_feats + 1)) + list(self.potential_labels))

            #print('computed label:', label)
            #print('self.potential_labels:', self.potential_labels)
            #print('self.potential_feats:', self.potential_feats)

            if self.options.mode == 'sep':
                self.selected_samps = self.samps_[label]
                self.non_selected_samps = {s for l in self.labels if l != label for s in self.samps_[l]}
            else:
                self.selected_samps = set(range(len(self.data.samps)))
                self.non_selected_samps = set()

            self.all_samps = self.selected_samps.union(self.non_selected_samps)

            wght = sum([self.data.wghts[i] for i in self.all_samps])

            if self.options.approx:
                self.lambda_ = int(math.ceil(wght * float(self.options.lambda_)))
            else:
                self.lambda_ = wght * decimal.Decimal(self.options.lambda_)

            if self.options.approach == 'sparse' and self.options.verb:
                print('c1 1 lit == {0} misses'.format(1 / self.lambda_))

            ubound = (self.nof_feats + 1) * wght

            # iterative over the number of literals
            self.nof_lits = 2

            while True:
                #print('label:', label)
                #print('self.nof_lits:', self.nof_lits)
                self.encode_constraints()

                with solver(self.enc, solver=self.options.solver,
                            adapt=self.options.am1, exhaust=self.options.exhaust,
                            minz=self.options.minz, trim=self.options.trim, ) as rc2:
                    model = rc2.compute()
                    if self.options.approach == 'sparse':
                        if model:
                            mis_samps = set(filter(lambda i: model[self.misclassified(i)-1] > 0, self.all_samps))
                            nof_miss = sum([self.data.wghts[i] for i in mis_samps])
                            nof_used = sum([1 if model[self.unused(j)-1] < 0 else 0 for j in range(1, self.nof_lits+1)])

                            if self.options.mode == 'sep':
                                self.selected_samps = self.samps_[label]
                                self.non_selected_samps = {s for l in self.labels if l != label for s in self.samps_[l]}
                            else:
                                self.selected_samps = set(range(len(self.data.samps)))
                                self.non_selected_samps = set()

                            # cost = nof_used + math.ceil(nof_miss / self.lambdas[label])
                            cost = nof_used + int(math.ceil(nof_miss / self.lambda_))
                            ubound = min(cost, ubound)

                            if self.options.verb:
                                print('c1 # of used: {0}; # of misses: {1} (out of {2}); ub: {3}'.format(nof_used, nof_miss, sum(self.data.wghts), ubound))

                            if ubound in (self.nof_lits, nof_used):
                                # either the formulation has reached the bound
                                # or the actual solution did
                                self.extract_cover(model)
                                self.nof_misses += nof_miss

                                for labl in self.samps_:
                                    if labl == label:
                                        self.samps_[labl] = set()
                                    else:
                                        self.samps_[labl] = self.samps_[labl].difference(mis_samps)
                                break
                            else:
                                if 10 < ubound - self.nof_lits:
                                    if 10 < 2 * self.nof_lits:
                                        self.nof_lits += 10
                                    else:
                                        self.nof_lits *= 2
                                else:
                                    self.nof_lits = ubound
                        else:
                            self.nof_lits += 1

                    else:
                        if model:
                            self.nof_misses = 0
                            self.extract_cover(model)
                            if self.options.mode == 'sep':
                                self.samps_[label] = set()
                            break
                        else:
                            self.nof_lits *= 2



        self.stime = resource.getrusage(resource.RUSAGE_SELF).ru_utime - self.init_stime
        self.ctime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime - self.init_ctime
        self.time = self.stime + self.ctime

        # computing accuracy
        nof_samps = sum(self.data.wghts)

        if not hasattr(self.data, 'wghts_filt'):
            self.accy_tot = 100. - (100. * self.nof_misses) / nof_samps
        else:
            self.accy = 100. - (100. * self.nof_misses) / nof_samps
            self.accy_tot = 100. - (100. * (self.nof_misses + self.data.wghts_filt)) / (
                        nof_samps + self.data.wghts_filt)

        #print('self.covrs:', self.covrs)

        if self.options.verb:
            size = len(self.covrs)
            wght = 0
            for rule in self.covrs:
                cond = []
                if rule['cond'] is not None:
                    for fv in rule['cond']:
                        feat = fv['feat']
                        value = fv['value']
                        sign = fv['sign']
                        cond.append("{0}'{1}: {2}'".format('not ' if not sign else '',
                                                           feat, value))
                else:
                    cond.append('true')
                wght += len(cond)
                label = rule['pred']['label']
                lv = rule['pred']['value']
                print("c1 cover: {0} => '{1}: {2}'".format(', '.join(cond), label, lv))

            if wght > 0:
                wght -= 1
            print('c2 cover nof_rules:', size)
            print('c2 cover nof_lits:', wght)
            print('c2 cover time: %.4f'%self.time)


        return self.covrs

    def unused(self, j):
        """
            True if literal at node j is unused.
        """

        return self.idpool.id('unused_{0}'.format(j))

    def feat(self, j, r):
        """
            True if literal at node j decides on feature r.
        """

        return self.idpool.id('feat_{0}_{1}'.format(j, r))

    def sign(self, j):
        """
            True if literal at node j is positive.
        """

        return self.idpool.id('sign_{0}'.format(j))

    def valid(self, i, j):
        """
            True iff sample i is valid at node j.
        """

        return self.idpool.id('valid_{0}_{1}'.format(i, j))

    def notclassified(self, i, j):
        """
            True if sample i is not previously classified by any node before j.
        """

        return self.idpool.id('notcls_{0}_{1}'.format(i, j))

    def misclassified(self, i):
        """
            True if sample i is misclassified.
        """

        return self.idpool.id('miscls_{0}'.format(i))

    def end(self, j):
        """
            True if node j is the end of the DL.
        """
        return self.idpool.id('end_{0}'.format(j))

    def aux0(self, j, r):
        """
            True iff literal at node j decides on feature r and it is negative.
        """
        return self.idpool.id('aux0_{0}_{1}'.format(j, r))

    def aux1(self, j, r):
        """
            True iff literal at node j decides on feature r and it is positive.
        """
        return self.idpool.id('aux1_{0}_{1}'.format(j, r))

    def auxvl(self, i, j):
        """
            True iff there exists node j is a leat node and sample i is valid at node j.
        """
        return self.idpool.id('auxvl_{0}_{1}'.format(i, j))

    def extract_cover(self, model):
        """
            Extracts a resulting DL from a model returned by a MaxSAT oracle.
        """
        label = self.data.names[-1]
        cond = []

        nodes = []
        unuse = ['{0}: {1}'.format(j, 'unused' if model[self.unused(j)-1] > 0 else 'used') for j in range(1, self.nof_lits + 1)]

        for j in range(1, self.nof_lits + 1):
            if model[self.unused(j)-1] > 0:
                break
            for r in self.potential_feats:
                if model[self.feat(j, r)-1] > 0:
                    nodes.append('node {0}: {1}'.format(j, r if model[self.sign(j)-1]>0 else -r))
                    if r > self.nof_feats:
                        # leaf node
                        if self.nof_labls >= 3:
                            value = self.data.fvmap.opp[r][-1]
                        else:
                            sign = model[self.sign(j)-1] > 0
                            if sign:
                                value = self.data.fvmap.opp[r][-1]
                            else:
                                value = self.data.fvmap.opp[r+1][-1]
                        pred = {'label': label, 'value': value}
                        cond_len = len(cond)
                        cond = cond if cond_len > 0 else None
                        self.covrs.append({'pred': pred, 'cond': cond})
                        self.cost += cond_len
                        cond = []
                    else:
                        feat, value = self.data.fvmap.opp[r]
                        if feat in value:
                            feat = value
                            value = 1
                        sign = model[self.sign(j)-1] > 0
                        cond.append({'feat': feat, 'value': value, 'sign': sign})
                    break

        #print(', '.join(unuse))
        #print('nodes:')
        #print(', '.join(nodes))
        return self.covrs

    def encode_constraints(self):
        """
            Encode the constraints
        """

        self.reset_idpool()

        for j in range(1, self.nof_lits + 1):
            for r in self.potential_feats:
                self.feat(j, r)
        for j in range(1, self.nof_lits + 1):
            self.sign(j)

        for j in range(1, self.nof_lits + 1):
            self.unused(j)

        self.enc = WCNF()

        # exactly one feature per node (including class/label)
        # Exactly one u[j] + S[j1] + S[j2]... + S[jf] + S[jc1] +  S[jc2] + ....
        #def cons_1(self):
        for j in range(1, self.nof_lits + 1):
            feats = [self.feat(j, r) for r in self.potential_feats] + [self.unused(j)]
            one = CardEnc.equals(lits=feats, vpool=self.idpool, encoding=self.options.enc)
            self.enc.extend(one)

        # Auxiliary: a0[j, r] <-> s[j, r] /\ -t[j], a1[j, r] <-> s[j, r] /\ t[j]
        #def cons_au(self):
        for j in range(1, self.nof_lits+1):
            for r in range(1, self.nof_feats+1):
                # Auxiliary: a0[j, r] <-> s[j, r] /\ -t[j]
                self.enc.append([-self.aux0(j, r), self.feat(j, r)])
                self.enc.append([-self.aux0(j, r), -self.sign(j)])
                self.enc.append([self.aux0(j, r), -self.feat(j, r), self.sign(j)])

                # Auxiliary: a1[j, r] <-> s[j, r] /\ t[j]
                self.enc.append([-self.aux1(j, r), self.feat(j, r)])
                self.enc.append([-self.aux1(j, r), self.sign(j)])
                self.enc.append([self.aux1(j, r), -self.feat(j, r), -self.sign(j)])


        # If a nodejis unused then so are all the following nodes
        # u[j] -> u[j+1]
        #def cons_2(self):
        self.enc.extend([[-self.unused(j), self.unused(j+1)] for j in range(1, self.nof_lits)])

        # The last used node is a leaf (cons_3, 4)
        # u[j+1] -> u[j]  \/  {for any c in C} s[j,c]
        #def cons_3(self):
        for j in range(1, self.nof_lits):
            f1 = [-self.unused(j+1), self.unused(j)] + \
                 [self.feat(j, c) for c in self.potential_labels]
            self.enc.append(f1)

        # u[N]  \/  {for any c in C} s[N,c]
        #def cons_4(self):
        j = self.nof_lits
        f1 = [self.unused(j)] + \
             [self.feat(j, c) for c in self.potential_labels]
        self.enc.append(f1)

        # all examples are not previously classified at the first node
        # n[i, 1] for all data
        #def cons_5(self):
        self.enc.extend([[self.notclassified(i, 1)] for i in self.all_samps])

        # All examples are valid at the first node
        # v[i,1]
        #def cons_7(self):
        self.enc.extend([[self.valid(i, 1)] for i in self.all_samps])

        # An example e_i is valid at nodej+ 1 iff j is a leaf node and it was previously unclassified,
        # or e_i is valid at node j and e_i and node j agree on the value of the feature s_jr selected for that node
        # v[i, j+1] <-> (s[j,c] /\ n[i,j+1]) \/ (v[i,j] /\ {for any r in K} (s[j,r] /\ (t[j] = pi_[i, r])))
        #def cons_8(self):
        #for i in range(1, self.nof_samps + 1):
        for i in self.all_samps:
            for j in range(1, self.nof_lits):
                self.enc.append([-self.valid(i, j+1), self.valid(i, j)] + \
                                 [self.feat(j, c) for c in self.potential_labels])
                self.enc.append([-self.valid(i, j+1), self.valid(i, j), self.notclassified(i, j+1)])
                f3 = []
                f3.append(-self.valid(i, j+1))
                for r in range(1, self.nof_feats+1):
                    if self.data.samps[i][r-1] < 0:
                        f3.append(self.aux0(j, r))
                    else:
                        f3.append(self.aux1(j, r))
                f3.extend([self.feat(j, c) for c in self.potential_labels])
                self.enc.append(f3)

                f4 = []
                f4.append(-self.valid(i, j+1))
                for r in range(1, self.nof_feats+1):
                    if self.data.samps[i][r-1] < 0:
                        f4.append(self.aux0(j, r))
                    else:
                        f4.append(self.aux1(j, r))
                f4.append(self.notclassified(i, j+1))
                self.enc.append(f4)

                disj = [self.valid(i, j+1), -self.notclassified(i, j+1)]
                self.enc.extend([disj + [-self.feat(j, c)]
                                 for c in self.potential_labels])

                disj = [self.valid(i, j+1), -self.valid(i, j)]
                for r in range(1, self.nof_feats+1):
                    f6 = disj[:]
                    if self.data.samps[i][r-1] < 0:
                        f6.append(-self.aux0(j, r))
                    else:
                        f6.append(-self.aux1(j, r))
                    self.enc.append(f6)

        # When there are 3 or more classes we restrict leaf nodes to only consider
        # true examples of the class
        # s[j,c] -> t[j]
        # def cons_10(self):
        if self.nof_labls >= 3:
            for j in range(1, self.nof_lits+1):
                self.enc.extend([[-self.feat(j, c), self.sign(j)]
                                 for c in self.potential_labels])

        elif self.options.mode == 'sep':
            # make sure the correct sign for separated perfect binary classification
            for j in range(1, self.nof_lits+1):
                if self.label == self.nof_feats + 1:
                    self.enc.append([-self.feat(j, self.nof_feats + 1), self.sign(j)])
                else:
                    self.enc.append([-self.feat(j, self.nof_feats + 1), -self.sign(j)])


        # An example e_i is previously unclassified at node j+ 1
        # iff it was previously unclassified, and either j is a not leaf node or it was invalid at the previous leaf node
        # (so not classified by the rule that finished there)
        # n[i, j+1] <->  n[i,j] /\ ((for and c in C, -s[j,c]) \/ -v[i,j])
        #def cons_6(self):
        def cons_classify():
            for i in self.all_samps:
                for j in range(1, self.nof_lits):
                    self.enc.append([-self.notclassified(i, j+1), self.notclassified(i, j)])
                    self.enc.extend([[-self.notclassified(i, j+1), -self.valid(i, j), -self.feat(j, c)]
                                     for c in self.potential_labels])
                    self.enc.append([self.notclassified(i,j+1), -self.notclassified(i,j)] +
                                    [self.feat(j, c) for c in self.potential_labels])
                    self.enc.append([self.notclassified(i,j+1), -self.notclassified(i,j), self.valid(i,j)])

        # self.M + self.non_M + 1)
        # If example e_i is valid at a leaf node j, they should agree on the class feature
        # s[j,c] /\ v[i, j] -> (t[j] = c[i])
        #def cons_9(self):
        def cons_valid_leaf():
            if self.nof_labls >= 3:
                for i in self.selected_samps:
                    for j in range(1, self.nof_lits+1):
                        for c in self.potential_labels:
                            if self.data.samps[i][-1] != c:
                                self.enc.append([-self.feat(j, c), -self.valid(i, j)])
            else:
                for i in self.selected_samps:
                    for j in range(1, self.nof_lits+1):
                        f = [-self.feat(j, self.nof_feats+1), -self.valid(i, j)]
                        if self.data.samps[i][-1] == self.nof_feats+1:
                            f.append(self.sign(j))
                        else:
                            f.append(-self.sign(j))
                        self.enc.append(f)


        # For every example there should be at least one leaf node where it is valid:
        # {For any j in N} ( {for any c in C} s[j,c] /\ v[i,j] )
        #def cons_11(self):
        def cons_cover():
            # Auxiliary variables vl[i,j] <-> {for any c in C }s[jc] /\ v[ij]
            for i in self.all_samps:
                for j in range(1, self.nof_lits+1):
                    self.enc.append([-self.auxvl(i, j)] +
                                    [self.feat(j, c) for c in self.potential_labels])
                    self.enc.append([-self.auxvl(i, j), self.valid(i, j)])
                    disj = [self.auxvl(i, j), -self.valid(i, j)]
                    self.enc.extend([disj + [-self.feat(j, c)]
                                     for c in self.potential_labels])

            # For selected items, at least valid at one leaf
            # {For all i in selected items} {For any j in N} (vl[i, j])
            for i in self.selected_samps:
                self.enc.append([self.auxvl(i, j) for j in range(1, self.nof_lits+1)])

            # For selected items, not valid at all leaves
            # For non-selected items, {For any j in N} ( {for any c in C} s[j,c] /\ v[i,j] )
            # - ({for any c in C} s[jc] /\ v[i, j])
            # (-vl[i,j])
            if self.options.mode == 'sep':
                for i in self.non_selected_samps:
                    self.enc.extend([[-self.auxvl(i, j)] for j in range(1, self.nof_lits+1)])

        # The last used node ends a decision list
        # u[j+1] /\ -u[j] -> x[j]
        #def cons_12(self):
        def cons_lastuse():
            for j in range(1, self.nof_lits):
                self.enc.append([-self.unused(j+1), self.unused(j), self.end(j)])

        # An end node is always a leaf
        # x[j] -> s[j,c1] \/ s[j,c2] ....
        #def cons_13(self):
        def cons_endleaf():
            for j in range(1, self.nof_lits+1):
                self.enc.append([-self.end(j)] +
                                [self.feat(j, c) for c in self.potential_labels])

        # An example e_i is previously unclassified at node j+ 1
        # iff j is an end of decision list node, or it was previously unclassified, and either j is a not leaf node
        # or it was invalid at the previous leaf node (so not classified by the rule that finished there)
        # n[i,j+1] <->x[j] \/ (n[i,j] /\( ({any c in C}-s[j,c]) \/ -v[i,j] ) )
        #def cons_14(self):
        def cons_classify_end():
            for i in self.all_samps:
                for j in range(1, self.nof_lits):
                    self.enc.append([-self.end(j), self.notclassified(i, j+1)])
                    self.enc.append([-self.notclassified(i, j), self.notclassified(i, j+1)] +
                                    [self.feat(j, c) for c in self.potential_labels])
                    self.enc.append([self.valid(i, j), -self.notclassified(i, j), self.notclassified(i, j+1)])
                    self.enc.append([self.notclassified(i, j), self.end(j), -self.notclassified(i, j+1)])
                    f5 = [-self.valid(i, j), self.end(j), -self.notclassified(i, j+1)]
                    self.enc.extend([f5 + [-self.feat(j, c)]
                                     for c in self.potential_labels])

        # If example e_i is valid at a leaf node j
        # then they agree on the class feature or the item is misclassified
        # s[j,c] /\ v[i, j] -> t[j] = c[i] \/ m[i]
        #def cons_15(self):
        def cons_misclassify():
            if self.nof_labls >= 3:
                # s[j,c] /\ v[i, j] -> c[i] \/ m[i]
                for i in self.selected_samps:
                    for j in range(1, self.nof_lits+1):
                        for c in self.potential_labels:
                            if self.data.samps[i][-1] != c:
                                self.enc.append([-self.feat(j, c), -self.valid(i, j), self.misclassified(i)])
            else:
                for i in self.selected_samps:
                    for j in range(1, self.nof_lits+1):
                        f = [-self.feat(j, self.nof_feats+1), -self.valid(i, j), self.misclassified(i)]
                        if self.data.samps[i][-1] == self.nof_feats+1:
                            f.append(self.sign(j))
                        else:
                            f.append(-self.sign(j))
                        self.enc.append(f)


        # For every example there should be at least one leaf node where it is valid:
        # m[i] \/ {For any j in N} ( {for any c in C} s[j,c] /\ v[i,j] )
        #def cons_16():
        def cons_classifiedonce():
            # Auxiliary variables vl[i,j] <-> {for any c in C }s[jc] /\ v[ij]
            for i in self.all_samps:
                for j in range(1, self.nof_lits+1):
                    self.enc.append([-self.auxvl(i, j)] +
                                    [self.feat(j, c)
                                     for c in self.potential_labels])
                    self.enc.append([-self.auxvl(i, j), self.valid(i, j)])

                    disj = [self.auxvl(i, j), -self.valid(i, j)]
                    self.enc.extend([disj + [-self.feat(j, c)]
                                     for c in self.potential_labels])

            # For selected items, at least valid at one leaf or misclassified
            # {For all i in selected items} m[i] \/ {For any j in N} (vl[i, j])
            for i in self.selected_samps:
                self.enc.append([self.misclassified(i)] +
                                [self.auxvl(i, j) for j in range(1, self.nof_lits+1)])

            # For selected items, not valid at all leaves or misclassified
            # For non-selected items, {For any j in N} ( {for any c in C} s[j,c] /\ v[i,j] )
            # mi \/ - ({for any c in C} s[jc] /\ v[i, j])
            # mi \/ (-vl[i,j])
            if self.options.mode == 'sep':
                for i in self.non_selected_samps:
                    for j in range(1, self.nof_lits+1):
                        self.enc.append([self.misclassified(i), -self.auxvl(i, j)])

        #####################

        # if it is a hybrid model
        if self.options.approach == 'hybrid' and not self.options.mode == 'sep':
            #cons_1(self)
            #cons_au(self)
            #cons_2(self)
            #cons_3(self)
            #cons_4(self)
            #cons_5(self)
            #cons_7(self)
            #cons_8(self)
            cons_valid_leaf()
            #cons_10(self)
            cons_cover()
            cons_lastuse()
            cons_endleaf()
            cons_classify_end()
            # Maximise u[j]
            for j in range(1, self.nof_lits+1):
                self.enc.append([self.unused(j)], weight=1)

        else:
            # if it is a perfect DL model
            if self.options.approach != 'sparse':
                #cons_1(self)
                #cons_au(self)
                #cons_2(self)
                #cons_3(self)
                #cons_4(self)
                #cons_5(self)
                cons_classify()
                #cons_7(self)
                #cons_8(self)
                if self.options.mode != 'sep':
                    #print('complete perfect DL')
                    cons_valid_leaf()
                else:
                    #print('sep perfect DL')
                    pass
                #cons_10(self)
                cons_cover()
                # Maximise u[j]
                for j in range(1, self.nof_lits + 1):
                    self.enc.append([self.unused(j)], weight=1)
            # if it is a sparse DL model
            else:
                #cons_1(self)
                #cons_au(self)
                #cons_2(self)
                #cons_3(self)
                #cons_4(self)
                #cons_5(self)
                cons_classify()
                #cons_7(self)
                #cons_8(self)
                #cons_10(self)
                if self.options.mode != 'sep':
                    #print('complete sparse DL')
                    cons_misclassify()
                else:
                    #print('separated sparse DL')
                    pass
                cons_classifiedonce()

                # objective function
                for i in self.all_samps:
                    self.enc.append([-self.misclassified(i)],
                                    weight=self.options.misw * self.data.wghts[i])

                for j in range(1, self.nof_lits + 1):
                    self.enc.append([self.unused(j)], weight=self.lambda_)
