#!/usr/bin/env python
#-*- coding:utf-8 -*-
##
## cover.py
##
##  Created on: Dec 5, 2017
##      Author: Alexey S. Ignatiev
##      E-mail: aignatiev@ciencias.ulisboa.pt
##

#
#==============================================================================
import collections
import math
import os
from pysat.card import *
from pysat.examples.lbx import LBX
from pysat.examples.rc2 import RC2
from pysat.formula import WCNFPlus
import resource
import socket
import six
from six.moves import range
import sys

# checking whether cplex is available
cplex_present = False
# cplex_present = True
# try:
#     import cplex
# except ImportError:
#     cplex_present = False

# checking whether gurobi is available
gurobi_present = True
try:
    import gurobipy as gurobi
except ImportError:
    gurobi_present = False


#
#==============================================================================
class Coverer(object):
    """
        IHS MaxSAT-based prime cover computation.
    """

    def __init__(self, data, primes, options):
        """
            Constructor.
        """

        self.init_stime = resource.getrusage(resource.RUSAGE_SELF).ru_utime
        self.init_ctime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime

        self.data = data
        self.primes = primes
        self.options = options

        # clusterizing original samples
        self.clmap = collections.defaultdict(lambda: [])
        for i, s in enumerate(self.data.samps):
            self.clmap[s[-1]].append(i)

            self.data.samps[i] = set(s[:-1])  # needed for the compute() function

        # casting primes to sets
        for out in six.iterkeys(self.primes):
            for i, p in enumerate(self.primes[out]):
                self.primes[out][i] = set(p)

    def compute(self):
        """
            Cover samples for all labels (separately).
        """

        self.covers = {}
        self.cost = 0
        print('c2 computing smallest prime cover')

        if self.options.cover == 'mxsat':
            self.compute_mxsat()
        elif self.options.cover == 'cplex':
            self.compute_cplex()
        elif self.options.cover in ('gurobi', 'grb'):
            self.compute_gurobi()

        # recording time
        self.stime = resource.getrusage(resource.RUSAGE_SELF).ru_utime - self.init_stime
        self.ctime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime - self.init_ctime
        self.time = self.stime + self.ctime
        return self.covers

    def compute_mxsat(self):
        """
            Cover samples for all labels using MaxSAT or MCS enumeration.
        """

        # traversing all sample clusters
        # (each cluster is hit/covered separately)
        for cid, cluster in six.iteritems(self.clmap):

            # skip not interesting classes
            # (decided by the command-line option '-c')
            if cid not in self.primes:
                continue

            # we model a set cover problem with MaxSAT
            formula = WCNFPlus()

            # hard part of the formula
            if self.options.accuracy == 100.0:
                for sid in cluster:  # for every sample in cluster
                    to_hit = []

                    # print('sample', self.data.samps[sid])
                    for pid, prime in enumerate(self.primes[cid]):
                        # print('prime', pid + 1, prime)
                        if prime.issubset(self.data.samps[sid]):
                            to_hit.append(pid + 1)

                    # print('to hit:', to_hit)
                    formula.append(to_hit)
            else:
                topv = len(self.primes[cid])
                allvars = []

                # hard clauses first
                for sid in cluster:  # for every sample in cluster
                    to_hit = []

                    for pid, prime in enumerate(self.primes[cid]):
                        if prime.issubset(self.data.samps[sid]):
                            to_hit.append(pid + 1)

                    topv += 1
                    allvars.append(topv)
                    formula.append([-topv] + to_hit)
                    for pid in to_hit:
                        formula.append([topv, -pid])

                # forcing at least the given percentage of samples to be covered
                cnum = int(math.ceil(self.options.accuracy * len(allvars) / 100.0))
                al = CardEnc.atleast(allvars, bound=cnum, top_id=topv, encoding=self.options.enc)
                if al:
                    for cl in al.clauses:
                        formula.append(cl)

            # soft clauses
            for pid in range(len(self.primes[cid])):
                formula.append([-pid - 1], weight=1)

            if self.options.weighted and not self.options.approx:
                # it is safe to add weights for all primes
                # because each prime covers at least one sample

                formula.wght = [len(prime) for prime in self.primes[cid]]

            if self.options.pdump:
                fname = 'cover{0}.{1}@{2}.wcnf'.format(cid, os.getpid(), socket.gethostname())
                formula.to_file(fname)

            # choosing the right solver
            if not self.options.approx:
                MaxSAT = RC2Stratified if self.options.blo else RC2
                hitman = MaxSAT(formula, solver=self.options.solver,
                        adapt=self.options.am1, exhaust=self.options.exhaust,
                        trim=self.options.trim, minz=self.options.minz)
            else:
                hitman = LBX(formula, use_cld=self.options.use_cld,
                        solver_name=self.options.solver, use_timer=False)

            # and the cover is...
            if not self.options.approx:
                cover = list(filter(lambda l: 0 < l <= len(self.primes[cid]) + 1, hitman.compute()))
                self.covers[cid] = cover

                self.cost += hitman.cost
            else:
                # approximating by computing a number of MCSes
                covers = []
                for i, cover in enumerate(hitman.enumerate()):
                    hitman.block(cover)
                    if self.options.weighted:
                        cost = sum([len(self.primes[cid][pid - 1]) for pid in cover])
                    else:
                        cost = len(cover)

                    covers.append([cover, cost])

                    if i + 1 == self.options.approx:
                        break

                self.covers[cid], cost = min(covers, key=lambda x: x[1])
                self.cost += cost

            if self.options.verb:
                for pid in self.covers[cid]:
                    premise = []
                    mapped_back = {}

                    for l in self.primes[cid][pid - 1]:
                        if l not in self.data.ovmap:
                            feat, val = self.data.fvmap.opp[abs(l)]
                            premise.append('{0}\'{1}: {2}\''.format('' if l > 0 else 'not ', feat, val))
                        else:
                            feat = self.data.ovmap[l][0]

                            if feat not in mapped_back:
                                mapped_back[feat] = set(self.data.ovmap[l][1:])
                            else:
                                mapped_back[feat].intersection_update(set(self.data.ovmap[l][1:]))
                                self.cost -= 1

                    for feat, vals in six.iteritems(mapped_back):
                        if len(vals) == 1:
                            premise.append('\'{0}: {1}\''.format(feat, vals.pop()))
                        elif vals:
                            premise.append('\'{0} in ({1})\''.format(feat, ', '.join([v for v in vals])))
                        else:  # this seems impossible
                            assert 0, 'empty value for feature {0}'.format(feat)

                    print('c2 cover: {0} => \'{1}: {2}\''.format(', '.join(premise), *self.data.fvmap.opp[cid]))

            hitman.delete()

    def compute_cplex(self):
        """
            Cover samples for all labels using CPLEX.
        """

        assert cplex_present, 'CPLEX is unavailable'

        # traversing all sample clusters
        # (each cluster is hit/covered separately)
        for cid, cluster in six.iteritems(self.clmap):

            # skip not interesting classes
            # (decided by the command-line option '-c')
            if cid not in self.primes:
                continue

            # initializing the solver
            hitman = cplex.Cplex()

            # turning logger off
            hitman.set_log_stream(None)
            hitman.set_error_stream(None)
            hitman.set_results_stream(None)
            hitman.set_warning_stream(None)

            # variables
            vnames = ['p_{0}'.format(i + 1) for i in range(len(self.primes[cid]))]
            prvars = hitman.variables.add(names=vnames, types='B' * len(vnames))

            # hard constraints
            for sid in cluster:  # for every sample in cluster
                to_hit = []

                # print('sample', self.data.samps[sid])
                for pid, prime in enumerate(self.primes[cid]):
                    # print('prime', pid + 1, prime)
                    if prime.issubset(self.data.samps[sid]):
                        to_hit.append(vnames[pid])

                hitman.linear_constraints.add(lin_expr=[[to_hit, [1] * len(to_hit)]], senses=['G'], rhs=[1], names=['sid{0}'.format(sid)])

            # optimization criterion
            if self.options.weighted:
                hitman.objective.set_linear([(vnames[pid], len(prime)) for pid, prime in enumerate(self.primes[cid])])
            else:
                hitman.objective.set_linear([(vnames[pid], 1) for pid, prime in enumerate(self.primes[cid])])

            hitman.objective.set_sense(hitman.objective.sense.minimize)

            if self.options.pdump:
                fname = 'cover{0}.{1}@{2}.lp'.format(cid, os.getpid(), socket.gethostname())
                hitman.write(fname)

            # and the cover is...
            hitman.solve()
            model, cover = hitman.solution, []
            for pid, v in enumerate(vnames, 1):
                if int(model.get_values(v)) > 0:
                    cover.append(pid)

            self.covers[cid] = cover
            self.cost += int(hitman.solution.get_objective_value())

            if self.options.verb:
                for pid in self.covers[cid]:
                    premise = []
                    mapped_back = {}

                    for l in self.primes[cid][pid - 1]:
                        if l not in self.data.ovmap:
                            feat, val = self.data.fvmap.opp[abs(l)]
                            premise.append('{0}\'{1}: {2}\''.format('' if l > 0 else 'not ', feat, val))
                        else:
                            feat = self.data.ovmap[l][0]

                            if feat not in mapped_back:
                                mapped_back[feat] = set(self.data.ovmap[l][1:])
                            else:
                                mapped_back[feat].intersection_update(set(self.data.ovmap[l][1:]))
                                self.cost -= 1

                    for feat, vals in six.iteritems(mapped_back):
                        if len(vals) == 1:
                            premise.append('\'{0}: {1}\''.format(feat, vals.pop()))
                        elif vals:
                            premise.append('\'{0} in ({1})\''.format(feat, ', '.join([v for v in vals])))
                        else:  # this seems impossible
                            assert 0, 'empty value for feature {0}'.format(feat)

                    print('c2 cover: {0} => \'{1}: {2}\''.format(', '.join(premise), *self.data.fvmap.opp[cid]))

            # cleaning up
            for sid in cluster:
                hitman.linear_constraints.delete('sid{0}'.format(sid))

    def compute_gurobi(self):
        """
            Cover samples for all labels using Gurobi.
        """

        assert gurobi_present, 'Gurobi is unavailable'

        # traversing all sample clusters
        # (each cluster is hit/covered separately)
        for cid, cluster in six.iteritems(self.clmap):

            # skip not interesting classes
            # (decided by the command-line option '-c')
            if cid not in self.primes:
                continue

            # a hack to disable license file logging
            stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')

            # initializing the solver
            hitman = gurobi.Model()

            # restoring sys.stdout
            sys.stdout = stdout

            # turning logger off
            hitman.Params.OutputFlag = 0
            hitman.Params.LogToConsole = 0

            # variables
            vnames = []
            prvars = []
            for i, prime in enumerate(self.primes[cid]):
                vnames.append('p_{0}'.format(i + 1))
                if self.options.weighted:
                    prvars.append(hitman.addVar(vtype=gurobi.GRB.BINARY,
                        name=vnames[i], obj=len(prime)))
                else:
                    prvars.append(hitman.addVar(vtype=gurobi.GRB.BINARY,
                        name=vnames[i], obj=1))

            # hard constraints
            for sid in cluster:  # for every sample in cluster
                to_hit = []

                # print('sample', self.data.samps[sid])
                for pid, prime in enumerate(self.primes[cid]):
                    # print('prime', pid + 1, prime)
                    if prime.issubset(self.data.samps[sid]):
                        to_hit.append(prvars[pid])

                hitman.addConstr(lhs=gurobi.quicksum(1 * v for v in to_hit),
                        sense=gurobi.GRB.GREATER_EQUAL,
                        rhs=1,
                        name='sid{0}'.format(sid))

            if self.options.pdump:
                fname = 'cover{0}.{1}@{2}.lp'.format(cid, os.getpid(), socket.gethostname())
                hitman.write(fname)

            # and the cover is...
            hitman.optimize()
            cover = []
            for pid, prime in enumerate(prvars):
                if int(prime.X) > 0:
                    cover.append(pid)

            self.covers[cid] = cover
            self.cost += int(hitman.objVal)

            if self.options.verb:
                for pid in self.covers[cid]:
                    premise = []
                    mapped_back = {}

                    for l in self.primes[cid][pid]:
                        if l not in self.data.ovmap:
                            feat, val = self.data.fvmap.opp[abs(l)]
                            premise.append('{0}\'{1}: {2}\''.format('' if l > 0 else 'not ', feat, val))
                        else:
                            feat = self.data.ovmap[l][0]

                            if feat not in mapped_back:
                                mapped_back[feat] = set(self.data.ovmap[l][1:])
                            else:
                                mapped_back[feat].intersection_update(set(self.data.ovmap[l][1:]))
                                self.cost -= 1

                    for feat, vals in six.iteritems(mapped_back):
                        if len(vals) == 1:
                            premise.append('\'{0}: {1}\''.format(feat, vals.pop()))
                        elif vals:
                            premise.append('\'{0} in ({1})\''.format(feat, ', '.join([v for v in vals])))
                        else:  # this seems impossible
                            assert 0, 'empty value for feature {0}'.format(feat)

                    print('c2 cover: {0} => \'{1}: {2}\''.format(', '.join(premise), *self.data.fvmap.opp[cid]))

            # cleaning up
            for sid in cluster:
                c = hitman.getConstrByName('sid{0}'.format(sid))
                hitman.remove(c)

            hitman.remove(hitman.getVars())
