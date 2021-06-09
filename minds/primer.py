#!/usr/bin/env python
#-*- coding:utf-8 -*-
##
## primer.py
##
##  Created on: Dec 4, 2017
##      Author: Alexey S. Ignatiev
##      E-mail: aignatiev@ciencias.ulisboa.pt
##

#
#==============================================================================
import collections
import decimal
from grbmaxsat import Gurobi
from lbxplus import LBXPlus
from mcslsplus import MCSlsPlus
import os
from pysat.card import *
from pysat.examples.rc2 import RC2, RC2Stratified
from pysat.formula import WCNFPlus
import resource
import socket
import six
from six.moves import range


#
#==============================================================================
class Primer(object):
    """
        MaxSAT/MCS-based prime implicant enumerator.
    """

    def __init__(self, data, options):
        """
            Constructor.
        """

        self.init_stime = resource.getrusage(resource.RUSAGE_SELF).ru_utime
        self.init_ctime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime

        # saving data
        self.data = data

        # binarizing the data properly
        for i in range(len(self.data.samps)):
            samp_bin, out = self.data.samps[i][:-1], self.data.samps[i][-1]
            for l in samp_bin:
                if l > 0:  # negative literal means that the feature is binary
                    name, lit = self.data.fvmap.opp[l]
                    j = self.data.nm2id[name]

                    if len(self.data.feats[j]) > 2:
                        samp_bin += [-self.data.fvmap.dir[(name, l)] for l in sorted(self.data.feats[j].difference(set([lit])))]

            self.data.samps[i] = samp_bin + [out]

        # saving options
        self.options = options

        # create a MaxSAT formula for prime enumeration
        self.prepare_formula()

        # create and initialize primer
        if self.options.indprime == False:
            self.init_solver()

    def prepare_formula(self):
        """
            Prepare a MaxSAT formula for prime enumeration.
        """

        # creating a formula
        self.formula = WCNFPlus()

        # formula's variables
        self.orig_vars = max(self.data.fvmap.opp.keys())
        self.formula.nv = self.orig_vars * 2

        # creating soft clauses and hard p-clauses
        # as well as a mapping between dual-rail variables and input variables
        self.drvmap = {}
        for v in range(1, self.orig_vars + 1):
            if v not in self.data.deleted:
                self.formula.soft.append([-v])
                self.formula.soft.append([-v - self.orig_vars])

                self.formula.hard.append([-v, -v - self.orig_vars])  # p clauses

                self.drvmap[v] = v
                self.drvmap[v + self.orig_vars] = -v

        self.formula.wght = [1 for cl in self.formula.soft]
        self.formula.topw = len(self.formula.soft) + 1

        # hard clauses
        for sample in self.data.samps:
            cl = list(map(lambda l: -l if l < 0 else l + self.orig_vars, sample[:-1])) + [sample[-1]]
            self.formula.hard.append(cl)

        # additional hard clauses encoding that
        # each prime implicant should be meaningful,
        # i.e. it should contain exactly one class label assigned
        outs = [self.data.fvmap.dir[(self.data.names[-1], l)] for l in sorted(self.data.feats[-1])]
        one_out = CardEnc().equals(lits=outs, top_id=2 * self.orig_vars, encoding=self.options.enc)

        # updating top id (take Tseitin variables into account)
        self.formula.nv = one_out.nv

        # copying cardinality constraint
        for cl in one_out.clauses:
            cl_dr = list(map(lambda x: -x + self.orig_vars if abs(x) <= self.orig_vars and x < 0 else x, cl))
            self.formula.hard.append(cl_dr)

        # add hard clauses filtering out unnecessary combinations
        if self.options.filter or self.options.approach == 'pgreedy':
            self.filter_noncovers()

        if self.options.pdump:
            fname = 'primes.{0}@{1}.wcnf'.format(os.getpid(), socket.gethostname())
            self.formula.to_file(fname)

        if self.options.verb:
            print('c1 formula: {0}v, {1}c ({2}h+{3}s)'.format(self.formula.nv,
                len(self.formula.hard) + len(self.formula.soft),
                len(self.formula.hard), len(self.formula.soft)))

    def init_solver(self):
        """
            Create an initialize a solver for prime enumeration.
        """

        # initializing prime enumerator
        if self.options.primer == 'lbx':
            self.mcsls = LBXPlus(self.formula, use_cld=self.options.use_cld,
                    solver_name=self.options.solver, get_model=True,
                    use_timer=False)
        elif self.options.primer == 'mcsls':
            self.mcsls = MCSlsPlus(self.formula, use_cld=self.options.use_cld,
                    solver_name=self.options.solver, get_model=True,
                    use_timer=False)
        elif self.options.primer in ('grb', 'gurobi'):
            self.rc2 = Gurobi(self.formula, self.options)
        else:  # sorted or maxsat
            MaxSAT = RC2Stratified if self.options.blo else RC2
            self.rc2 = MaxSAT(self.formula, solver=self.options.solver,
                    adapt=self.options.am1, exhaust=self.options.exhaust,
                    trim=self.options.trim, minz=self.options.minz)

    def filter_noncovers(self):
        """
            Add hard constraints to block all non-covering primes.
        """

        topv = self.formula.nv
        ncls = len(self.formula.hard)
        tvars = []  # auxiliary variables

        allv = []
        for v in range(1, self.data.fvars + 1):
            allv.append(v)
            allv.append(v + self.orig_vars)
        allv = set(allv)

        self.tcls = {self.data.fvmap.dir[(self.data.names[-1], c)]: [] for c in sorted(self.data.feats[-1])}

        # a class to compute
        labels = None  # compute all
        if self.options.to_compute not in ('best', 'all'):
            to_compute = self.options.to_compute.split(',')
            labels = [self.data.fvmap.dir[self.data.names[-1], c] for c in to_compute]

        for sample in self.data.samps:
            if labels and sample[-1] not in labels:
                continue  # we are not interested in computing this class

            s = set([l if l > 0 else -l + self.orig_vars for l in sample[:-1]])

            # computing the complement of the sample
            compl = allv.difference(s)

            # encoding the complement (as a term) into a set of clauses
            if compl:
                topv += 1
                tvars.append(topv)
                self.tcls[sample[-1]].append(topv)

                compl = sorted(compl)
                for l in compl:
                    self.formula.hard.append([-l, -topv])

                self.formula.hard.append(compl + [topv])

        if tvars:
            self.tvars = tvars
            if self.options.approach == 'pgreedy':
                # replacing soft clauses for the greedy approach
                self.formula.soft = [[t] for t in self.tvars]
                self.formula.wght = [1 for t in self.tvars]
                self.formula.topw = len(self.formula.soft) + 1

            if not self.options.indprime:
                # add final clause forcing to cover at least one sample
                self.formula.hard.append(self.tvars[:])

                if self.options.plimit:
                    self.nof_p = {t: 0 for t in self.tvars}

            if self.options.verb:
                print('c1 blocked all non-covering primes')
                print('c1 added more {0} vars and {1} clauses'.format(
                    topv - self.formula.nv, len(self.formula.hard) - ncls))

            self.formula.nv = topv

    def compute(self):
        """
            Enumerate all prime implicants.
        """

        if self.options.primer in ('lbx', 'mcsls'):
            return self.compute_mcsls()
        else:  # sorted or maxsat
            return self.compute_sorted()

    def compute_mcsls(self):
        """
            Call an MCS enumerator.
        """

        print('c1 enumerating primes (mcs-based)')

        self.primes = collections.defaultdict(lambda: [])

        if self.options.approach == 'pgreedy':
            to_cover = set(self.tvars)

            while to_cover:
                covered, i = set(), 0
                for i, mcs in enumerate(self.mcsls.enumerate()):
                    mod = self.mcsls.get_model()
                    mcs = list(filter(lambda l: l > 0 and abs(l) <= 2 * self.orig_vars, mod))

                    prime, out = self.process_mcs(mcs)

                    # recording prime implicant
                    self.primes[out].append(prime)

                    # block
                    self.mcsls.add_clause([-l for l in mcs])

                    # samples that are covered by this last prime/mcs
                    prime_covered = set(self.tvars).intersection(set(mod))

                    # recording samples that are already covered
                    covered = covered.union(prime_covered)

                    if self.options.bsymm:
                        # breaking symmetric solutions
                        self.mcsls.add_clause(sorted(set(self.tvars).difference(prime_covered)))
                    elif i >= self.options.plimit and self.options.plimit > 0:
                        # filtering out already covered samples if necessary
                        to_cover = to_cover.difference(covered)
                        self.mcsls.add_clause(sorted(to_cover))
                        i, move_on = 0, False

                        if self.options.verb > 1:
                            print('c1 uncovered:', len(to_cover))

                if not i:
                    break

            self.mcsls.delete()
        elif self.options.indprime:
            self.formula.hard.append([])

            mcses = []
            for t in self.tvars:
                # set sample to cover
                self.formula.hard[-1] = [t]

                # initialize solver from scratch
                self.init_solver()

                # block all previosly computed MCSes
                for mcs in mcses:
                    self.mcsls.block(mcs)

                # enumerate MCSes covering sample t
                for i, mcs in enumerate(self.mcsls.enumerate()):
                    mod = self.mcsls.get_model()  # mcs will be extracted from a model
                    mcs = list(filter(lambda l: l > 0 and abs(l) <= 2 * self.orig_vars, mod))
                    prime, out = self.process_mcs(mcs)

                    # recording prime implicant
                    self.primes[out].append(prime)

                    # block
                    self.mcsls.add_clause([-l for l in mcs])

                    # record mcs to block later
                    mcses.append(mcs)

                    if self.options.bsymm:
                        # breaking symmetric solutions
                        self.mcsls.add_clause(sorted(set(self.tvars).difference(set(mod))))

                    if self.options.plimit and i + 1 == self.options.plimit:
                        break

                self.mcsls.delete()
        else:
            for mcs in self.mcsls.enumerate():
                mod = self.mcsls.get_model()
                mcs = list(filter(lambda l: l > 0 and abs(l) <= 2 * self.orig_vars, mod))

                prime, out = self.process_mcs(mcs)

                # recording prime implicant
                self.primes[out].append(prime)

                # block
                self.mcsls.add_clause([-l for l in mcs])

                if self.options.bsymm:
                    # breaking symmetric solutions
                    self.mcsls.add_clause(sorted(set(self.tvars).difference(set(mod))))

                # check if there are enough MCSes
                if self.options.plimit:
                    model = self.mcsls.get_model()

                    i, reduced = 0, False
                    while i < len(self.tvars):
                        t = self.tvars[i]
                        if model[t - 1] > 0:
                            self.nof_p[t] += 1

                        if self.nof_p[t] < self.options.plimit:
                            i += 1
                        else:
                            self.tvars[i] = self.tvars[-1]
                            self.tvars.pop()
                            reduced = True

                    if reduced:
                        self.mcsls.oracle.add_clause(self.tvars)

                        if not self.tvars:
                            break

            self.mcsls.delete()

        # recording time
        self.stime = resource.getrusage(resource.RUSAGE_SELF).ru_utime - self.init_stime
        self.ctime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime - self.init_ctime
        self.time = self.stime + self.ctime

        return self.primes

    def compute_sorted(self):
        """
            MaxSAT-based prime implicant enumeration.
        """

        print('c1 enumerating primes (maxsat-based)')

        self.primes = collections.defaultdict(lambda: [])

        self.mcses = []

        if self.options.approach == 'pgreedy':
            to_cover = set(self.tvars)

            while to_cover:
                covered = set()
                best, move_on, i = -1, False, 0
                for mod in self.rc2.enumerate():
                    mcs = list(filter(lambda l: l > 0 and abs(l) <= 2 * self.orig_vars, mod))
                    i += 1

                    # decide whether we need to focus deeper on other samples
                    cost = self.rc2.cost
                    if self.options.plimit == -1 and cost > best:
                        if best > 0:
                            move_on = True

                        best = cost

                    prime, out = self.process_mcs(mcs)

                    # recording prime implicant
                    self.primes[out].append(prime)

                    # block
                    self.rc2.add_clause([-l for l in mcs])

                    # recording the mcs for future blocking
                    self.mcses.append(mcs)

                    # samples that are covered by this last prime/mcs
                    prime_covered = set(self.tvars).intersection(set(mod))

                    # recording samples that are already covered
                    covered = covered.union(prime_covered)

                    if self.options.bsymm:
                        # breaking symmetric solutions
                        self.rc2.add_clause(sorted(set(self.tvars).difference(prime_covered)))
                    elif move_on or i >= self.options.plimit:
                        # filtering out already covered samples if necessary
                        to_cover = to_cover.difference(covered)
                        self.rc2.add_clause(sorted(to_cover))
                        i, move_on = 0, False

                        if self.options.verb > 1:
                            print('c1 uncovered:', len(to_cover))

                if not i:
                    break

            self.rc2.delete()
        elif self.options.indprime:
            self.formula.hard.append([])

            mcses = []
            for t in self.tvars:
                # set sample to cover
                self.formula.hard[-1] = [t]

                # initialize solver from scratch
                self.init_solver()

                # block all previosly computed MCSes
                for mcs in mcses:
                    self.rc2.add_clause([-l for l in mcs])

                # enumerate MCSes covering sample t
                for i, mod in enumerate(self.rc2.enumerate()):
                    mcs = list(filter(lambda l: l > 0 and abs(l) <= 2 * self.orig_vars, mod))

                    # blocking the mcs properly
                    self.rc2.add_clause([-l for l in mcs])

                    # processing it
                    prime, out = self.process_mcs(mcs)

                    # recording prime implicant
                    self.primes[out].append(prime)

                    if self.options.bsymm:
                        # breaking symmetric solutions
                        self.rc2.add_clause(sorted(set(self.tvars).difference(set(mod))))

                    # record mcs to block later
                    mcses.append(mcs)

                    if self.options.plimit and i + 1 == self.options.plimit:
                        break

                self.rc2.delete()
        else:
            for mod in self.rc2.enumerate():
                mcs = list(filter(lambda l: l > 0 and abs(l) <= 2 * self.orig_vars, mod))

                # blocking the mcs properly
                self.rc2.add_clause([-l for l in mcs])

                # processing it
                prime, out = self.process_mcs(mcs)

                # recording the mcs for future blocking
                self.mcses.append(mcs)

                # recording prime implicant
                self.primes[out].append(prime)

                if self.options.bsymm:
                    # breaking symmetric solutions
                    self.rc2.add_clause(sorted(set(self.tvars).difference(set(mod))))

                # check if there are enough MCSes
                if self.options.plimit:
                    model = self.rc2.model

                    i, reduced = 0, False
                    while i < len(self.tvars):
                        t = self.tvars[i]
                        if model[t - 1] > 0:
                            self.nof_p[t] += 1

                        if self.nof_p[t] < self.options.plimit:
                            i += 1
                        else:
                            self.tvars[i] = self.tvars[-1]
                            self.tvars.pop()
                            reduced = True

                    if reduced:
                        self.rc2.add_clause(self.tvars)

                        if not self.tvars:
                            break

            self.rc2.delete()

        # recording time
        self.stime = resource.getrusage(resource.RUSAGE_SELF).ru_utime - self.init_stime
        self.ctime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime - self.init_ctime
        self.time = self.stime + self.ctime

        return self.primes

    def process_mcs(self, mcs):
        """
            Extract a prime implicant from MCS.
        """

        prime, out = [], None

        for i in mcs:
            # getting the corresponding variable
            v_orig = self.drvmap[i]

            # filtering out all labels
            if self.data.nm2id[self.data.fvmap.opp[abs(v_orig)][0]] == len(self.data.names) - 1:
                if v_orig > 0:
                    out = v_orig

                continue

            prime.append(v_orig)

        # printing prime implicant
        if self.options.verb > 1:
            premise = []

            for l in prime:
                feat, val = self.data.fvmap.opp[abs(l)]
                premise.append('{0}\'{1}: {2}\''.format('' if l > 0 else 'not ', feat, val))

            if self.options.verb > 2:
                print('c1 mcs: {0}'.format(' '.join([str(l) for l in mcs])))

            print('c1 prime: {0} => \'{1}: {2}\''.format(', '.join(premise), *self.data.fvmap.opp[out]))

        return prime, out


#
#==============================================================================
class PrimerCG(Primer, object):
    """
        MaxSAT/MCS-based prime implicant enumerator used in column generation.
    """

    def __init__(self, data, dobj, cid, primer, computed, options, node=None,
            left=False):
        """
            Constructor.
        """

        self.init_stime = resource.getrusage(resource.RUSAGE_SELF).ru_utime
        self.init_ctime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime

        # copying the parameters from the original primer
        self.data = data
        self.tcls = primer.tcls
        self.labl = cid
        self.options = options
        self.drvmap = primer.drvmap
        self.formula = primer.formula.copy()
        self.orig_vars = primer.orig_vars

        # modifying the soft clauses
        self.formula.soft = []
        self.formula.wght = []

        # dual objective coefficients
        self.dobj = dobj[:]

        # the positive part of the objective
        for v in range(1, self.orig_vars + 1):
            if v not in self.data.deleted:
                self.formula.soft.append([-v])
                self.formula.soft.append([-v - self.orig_vars])

        self.formula.wght = [1 for cl in self.formula.soft]

        # adding the negative part
        for i, (d, t) in enumerate(zip(self.dobj, self.tcls[self.labl])):
            if d:
                self.formula.soft.append([t])
                self.formula.wght.append(d)

        self.formula.topw = sum(self.formula.wght) + 1

        if node:
            t1, t2 = self.tcls[self.labl][node[0]], self.tcls[self.labl][node[1]]

            if left:
                # t1 and t2 must be equivalent, i.e.
                # the rows must be covered or uncovered together
                self.formula.append([-t1, t2])
                self.formula.append([t1, -t2])
            else:
                # t1 and t2 must be opposite
                self.formula.append([-t1, -t2])
                self.formula.append([ t1,  t2])

        if self.options.pdump:
            fname = 'primes-cg.{0}@{1}.wcnf'.format(os.getpid(), socket.gethostname())
            self.formula.to_file(fname)

        # already computed primes
        self.computed = computed

        # make sure we use an exact MaxSAT solver
        if self.options.primer in ('lbx', 'mcsls'):
            self.options.primer = 'sorted'

        # initializing the solver with the constructed formula
        self.init_solver()

        # blocking all previously found solutions
        for mcs in self.computed:
            self.rc2.add_clause([-l for l in mcs])

        # forcing the right class
        self.rc2.add_clause([self.labl])

        # list of primes
        self.primes = collections.defaultdict(lambda: [])

    def compute(self):
        """
            Compute a new prime.
        """

        model = self.rc2.compute()

        # stop if there is no other MaxSAT solution
        if not model:
            return

        # extracting the corresponding MCS
        mcs = list(filter(lambda l: l > 0 and abs(l) <= 2 * self.orig_vars, model))

        # blocking the mcs properly
        self.rc2.add_clause([-l for l in mcs])

        # recording the mcs for future blocking
        self.computed.append(mcs)

        # checking whether the positive or
        # the negative part of the cost is greater
        if self.options.weighted:
            # positive part (we need to filter out all label literals)
            pos = len(list(filter(lambda i: self.data.nm2id[self.data.fvmap.opp[abs(self.drvmap[i])][0]] < len(self.data.names) - 1, mcs)))
        else:
            # if the problem is unweighted, every prime has weight 1
            pos = 1

        # negative part (counting the negative cost)
        neg = decimal.Decimal(0) if self.options.primer == 'sorted' else 0
        for i, t in enumerate(self.tcls[self.labl]):
            if model[t - 1] > 0:
                neg += self.dobj[i]

        # the prime is interesting only if the total cost is negative
        if pos < neg:
            # processing it
            prime, out = self.process_mcs(mcs)

            # recording prime implicant
            self.primes[out].append(prime)

            return prime

    def process_mcs(self, mcs):
        """
            Extract a prime implicant from MCS.
        """

        prime, out = [], None

        for i in mcs:
            # getting the corresponding variable
            v_orig = self.drvmap[i]

            # filtering out all labels
            if self.data.nm2id[self.data.fvmap.opp[abs(v_orig)][0]] == len(self.data.names) - 1:
                if v_orig > 0:
                    out = v_orig

                continue

            prime.append(v_orig)

        # printing prime implicant
        if self.options.verb > 1:
            premise = []

            for l in prime:
                feat, val = self.data.fvmap.opp[abs(l)]
                premise.append('{0}\'{1}: {2}\''.format('' if l > 0 else 'not ', feat, val))

            if self.options.verb > 2:
                print('c2 mcs: {0}'.format(' '.join([str(l) for l in mcs])))

            print('c2 prime: {0} => \'{1}: {2}\''.format(', '.join(premise), *self.data.fvmap.opp[out]))

        return prime, out
