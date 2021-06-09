#!/usr/bin/env python
#-*- coding:utf-8 -*-
##
## hcover.py
##
##  Created on: Apr 6, 2020
##      Author: Alexey Ignatiev
##      E-mail: alexey.ignatiev@monash.edu
##

#
#==============================================================================
import collections
import decimal
from enum import Enum
import functools
import itertools
import math
import numpy as np
import os
from primer import PrimerCG
import resource
import socket
import six
from six.moves import range
import sys

# checking whether gurobi is available
gurobi_present = True
try:
    import gurobipy as gurobi
except ImportError:
    gurobi_present = False


#
#==============================================================================
class Status(Enum):
    """
        Status of a branch.
    """

    INFEASIBLE = 0
    FRACTIONAL = 1
    INTEGRAL = 2


#
#==============================================================================
class HybCoverer(object):
    """
        Hybrid coverer based on column generation.
    """

    def __init__(self, data, primer, primes, options):
        """
            Constructor.
        """

        assert gurobi_present, 'Gurobi is unavailable'

        self.init_stime = resource.getrusage(resource.RUSAGE_SELF).ru_utime
        self.init_ctime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime

        self.data = data
        self.primes = primes
        self.options = options

        # original primer used to generate an over-approximation
        self.prorig = primer

        # clusterizing original samples
        self.clmap = collections.defaultdict(lambda: [])
        for i, s in enumerate(self.data.samps):
            self.clmap[s[-1]].append(i)

            self.data.samps[i] = set(s[:-1])  # needed for the compute() function

        # casting primes to sets
        for out in six.iterkeys(self.primes):
            for i, p in enumerate(self.primes[out]):
                self.primes[out][i] = set(p)

        # initially the "set cover" solver is none
        self.hitman = None
        self.constr = []
        self.dtrail = []
        self.ltrail = []

        self.cost = 0

    def compute(self):
        """
            Cover samples for all labels (separately).
        """

        self.covers = {}
        self.costs = {}
        print('c2 computing smallest prime cover (based on column generation)')

        # covering every class independently
        for cid, cluster in six.iteritems(self.clmap):
            # skip not interesting classes
            # (decided by the command-line option '-c')
            if cid not in self.primes:
                continue

            # cover the target class
            self.cover_class(cid, cluster)

        self.cost = sum(self.costs.values())

        # recording time
        self.stime = resource.getrusage(resource.RUSAGE_SELF).ru_utime - self.init_stime
        self.ctime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime - self.init_ctime
        self.time = self.stime + self.ctime

        return self.covers

    def cover_class(self, cid, cluster):
        """
            Do branch-and-price for computing a cover for this class.
        """

        # previously computed MCSes
        self.mcses = self.prorig.mcses[:]

        # initialize the set cover solver
        self.init_hitman()

        # initial cover
        self.initial_cover(cid, cluster)

        # first, considering root node
        if self.node_cover(cid, cluster) == Status.FRACTIONAL:
            # if no integer solution at the root node, we need to branch
            print('c2 no integral solution at root node')
            print('c2 do branching')

            self.next_branch(cid, cluster)

        # report solution
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
                            self.costs[cid] -= 1

                for feat, vals in six.iteritems(mapped_back):
                    if len(vals) == 1:
                        premise.append('\'{0}: {1}\''.format(feat, vals.pop()))
                    elif vals:
                        premise.append('\'{0} in ({1})\''.format(feat, ', '.join([v for v in vals])))
                    else:  # this seems impossible
                        assert 0, 'empty value for feature {0}'.format(feat)

                print('c2 cover: {0} => \'{1}: {2}\''.format(', '.join(premise), *self.data.fvmap.opp[cid]))

    def init_hitman(self):
        """
            Initialize a new set cover solver.
        """

        # cleaning up
        if self.hitman:
            for c in self.constr:
                self.hitman.remove(c)

            self.constr = []
            self.hitman.remove(self.hitman.getVars())
            self.dobj = []
        else:
            # a hack to disable license file logging
            stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')

            # initializing the solver
            self.hitman = gurobi.Model()

            # restoring sys.stdout
            sys.stdout = stdout

            # turning logger off
            self.hitman.Params.OutputFlag = 0
            self.hitman.Params.LogToConsole = 0

        # empty cover matrix
        self.matrix = []

        # already considered pairs to branch on
        self.bpairs = set()

    def initial_cover(self, cid, cluster):
        """
            Compute an initial cover for the given class.
        """

        # variables
        self.vnames = []
        self.prvars = []
        for i, prime in enumerate(self.primes[cid]):
            self.vnames.append('p_{0}'.format(i + 1))

            if self.options.weighted:
                self.prvars.append(self.hitman.addVar(name=self.vnames[i], obj=len(prime)))
            else:
                self.prvars.append(self.hitman.addVar(name=self.vnames[i], obj=1))

        # hard constraints
        for sid in cluster:  # for every sample in cluster
            to_hit, row = [], []

            for pid, prime in enumerate(self.primes[cid]):
                if prime.issubset(self.data.samps[sid]):
                    to_hit.append(self.prvars[pid])
                    row.append(True)
                else:
                    row.append(False)

            c = self.hitman.addConstr(functools.reduce(lambda x, y: x + y, to_hit) >= 1)
            self.constr.append(c)
            self.matrix.append(row)

        # casting the matrix to be a numpy array
        self.matrix = np.array(self.matrix, dtype=bool)

        if self.options.pdump:
            fname = 'cover{0}-init.{1}@{2}.lp'.format(cid, os.getpid(), socket.gethostname())
            self.hitman.write(fname)

        # and the initial cover is...
        self.hitman.optimize()

        # dual objective coefficients
        self.dobj = []

        if self.options.primer == 'sorted':
            for c in self.constr:
                self.dobj.append(decimal.Decimal(c.Pi))
        else:
            for c in self.constr:
                self.dobj.append(c.Pi)

        if self.options.verb >= 3:
            print('c2 dual coeff:', ', '.join([str(v) for v in self.dobj]))

    def node_cover(self, cid, cluster, node=None, left=False):
        """
            Do column generation at the root node of the search tree.
        """

        self.primer = PrimerCG(self.data, self.dobj, cid, self.prorig,
                self.mcses, self.options, node=node, left=left)

        # generating primes incrementally until we get a negative cost
        nof_primes, i = len(self.primes[cid]), 0
        while True:
            prime = self.primer.compute()

            if not prime:
                break

            # add a new column
            self.append_column(set(prime), cluster)

            # record the new prime
            self.primes[cid].append(set(prime))

            i += 1
            if self.options.update and i % self.options.update == 0:
                # updating the dual objective coefficients
                self.hitman.optimize()
                if self.options.primer == 'sorted':
                    dobj_new = [decimal.Decimal(c.Pi) for c in self.constr]
                else:
                    dobj_new = [c.Pi for c in self.constr]

                if self.dobj != dobj_new:
                    self.dobj = dobj_new

                    # restarting the primer
                    self.primer = PrimerCG(self.data, self.dobj, cid, self.prorig,
                            self.mcses, self.options, node=node,
                            left=left)

                if self.options.verb >= 3:
                    print('c2 dual coeff:', ', '.join([str(v) for v in self.dobj]))

        print('c2 # of additional primes: {0}'.format(len(self.primes[cid]) - nof_primes))

        # trying to get an optimal solution
        self.hitman.optimize()

        if self.options.pdump:
            fname = 'cover{0}-i{3}{4}.{1}@{2}.lp'.format(cid, os.getpid(), socket.gethostname(), len(self.primes[cid]), left)
            self.hitman.write(fname)

        if self.options.verb >= 3:
            print('c2 hitman status:', self.hitman.status)

        if self.hitman.status == gurobi.GRB.INFEASIBLE:
            return Status.INFEASIBLE

        if self.options.verb >= 3:
            print('c2 objective:', self.hitman.objVal)

        cover = []
        is_integral = True
        for pid, prime in enumerate(self.prvars):
            if not prime.X.is_integer():
                is_integral = False
                break

            if int(prime.X) > 0:
                cover.append(pid)

        if is_integral:
            # integral solution is obtained => save it for later
            if not cid in self.costs or self.hitman.objVal < self.costs[cid]:
                self.covers[cid] = cover
                self.costs[cid] = int(self.hitman.objVal)

                if self.options.verb >= 3:
                    print('c2 improved solution found')

            return Status.INTEGRAL
        else:
            if cid in self.costs and self.hitman.objVal >= self.costs[cid]:
                # technically, we should return FRACTIONAL but I am lazy...
                return Status.INFEASIBLE

            self.nonint = []

            for pid, prime in enumerate(self.prvars):
                if not prime.X.is_integer():
                    self.nonint.append(pid)

            return Status.FRACTIONAL

    def append_column(self, prime, cluster):
        """
            Add a new column on demand.
        """

        # characteristic function
        indicators = []

        for sid in cluster:  # for every sample in cluster
            indicators.append(1 if prime.issubset(self.data.samps[sid]) else 0)

        column = gurobi.Column(indicators, self.constr)

        # adding the column to the matrix
        self.matrix = np.append(self.matrix,
                np.array([indicators], dtype=bool).T, axis=1)

        # name of the new variable
        pid = len(self.vnames)
        self.vnames.append('p_{0}'.format(pid + 1))

        if self.options.weighted:
            self.prvars.append(self.hitman.addVar(name=self.vnames[pid],
                obj=len(prime), column=column))
        else:
            self.prvars.append(self.hitman.addVar(name=self.vnames[pid],
                obj=1, column=column))

    def branch_cover(self, cid, cluster, node, left=True):
        """
            Do Ryan-Foster branching.
        """

        if self.options.verb >= 4:
            print('c2 matrix:')
            print(self.matrix)

        if self.options.verb >= 3:
            print('c2 branching on: {0} ({1})'.format(node, 'left' if left else 'right'))

        # new trail limit
        self.ltrail.append(len(self.dtrail))

        # column filtering based on the pair we currently branch on
        filt = lambda x, y: x != y if left else x == y
        for pid, prime in enumerate(self.prvars):
            if filt(*list(map(lambda r: self.matrix[r, pid], node))):
                self.dtrail.append(self.hitman.addConstr(prime == 0.0, '{0}_deleted'.format(self.vnames[pid])))

                if self.options.verb >= 3:
                    print('c2 adding constraints:', '{0}_deleted'.format(self.vnames[pid]))

        # here we branch further if a fractinal solution is
        # obtained & it is better than the current optimum
        if self.node_cover(cid, cluster, node=node, left=left) == Status.FRACTIONAL:
            self.next_branch(cid, cluster)

        # returning the filtered columns when backtracking
        while len(self.dtrail) != self.ltrail[-1]:
            if self.options.verb >= 3:
                print('c2 removing constraints:', self.dtrail[-1].getAttr('constrName'))

            self.hitman.remove(self.dtrail.pop())

        self.ltrail.pop()

    def next_branch(self, cid, cluster):
        """
            Compute the next node to branch on.
        """

        wmtx, branch_on = self.compute_wmatrix(), None
        for item in sorted(wmtx.items(), key=lambda item: item[1]):
            if item[0] in self.bpairs:
                continue

            branch_on = item[0]
            break

        if branch_on:
            self.bpairs.add(branch_on)

            # do recursive branching
            self.branch_cover(cid, cluster, branch_on, left=True)
            self.branch_cover(cid, cluster, branch_on, left=False)

    def compute_wmatrix(self):
        """
            Compute the weight matrix given the last solution.
        """

        wght = collections.defaultdict(lambda: 0.)

        # iterating over non-integer columns
        for col in self.nonint:
            val = self.prvars[col].X  # this should be a non-integer value

            # filter out zero rows in the current column
            rows = list(filter(lambda r: self.matrix[r, col], range(self.matrix.shape[0])))

            for i, row1 in enumerate(rows):
                for j in range(i + 1, len(rows)):
                    wght[(row1, rows[j])] += val

        return wght
