#!/usr/bin/env python
#-*- coding:utf-8 -*-
##
## grbmaxsat.py
##
##  Created on: Apr 14, 2020
##      Author: Alexey Ignatiev
##      E-mail: alexey.ignatiev@monash.edu
##

#
#==============================================================================
from __future__ import print_function
import functools
import os
from pysat.formula import WCNF
from six.moves import range
import socket
import sys

# checking whether gurobi is available
gurobi_present = True
try:
    import gurobipy as gp
    from gurobipy import GRB
except ImportError:
    gurobi_present = False


#
#==============================================================================
class Gurobi(object):
    """
        Gurobi-based oracle for MaxSAT replicating the interface of RC2.
    """

    def __init__(self, formula, options):
        """
            Constructor.
        """

        assert gurobi_present, 'Gurobi is unavailable'

        self.options = options

        # a hack to disable license file logging
        stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

        # initializing the solver
        self.oracle = gp.Model()

        # restoring sys.stdout
        sys.stdout = stdout

        # turning logger off
        self.oracle.Params.OutputFlag = 0
        self.oracle.Params.LogToConsole = 0

        # setting the sense to maximization
        self.oracle.setAttr('ModelSense', -1)

        # a mapping from the original MaxSAT
        # variables to internal Gurobi variables
        self.vmap = {}

        # constraints
        self.constr = []

        # feeding the oracle with the input MaxSAT formula
        self.add_formula(formula)

        if self.options.pdump:
            fname = 'primes.{0}@{1}.lp'.format(os.getpid(), socket.gethostname())
            self.oracle.write(fname)

    def add_formula(self, formula):
        """
            Add a formula to the solver.
        """

        # objective function
        for cl, w in zip(formula.soft, formula.wght):
            lit = cl[0]
            wght = w if lit > 0 else -w

            v = self.oracle.addVar(vtype=GRB.BINARY,
                    name='x{0}'.format(abs(lit)), obj=wght)

            self.vmap[abs(lit)] = v

        obj = self.oracle.getObjective()

        # hard clauses
        for cl in formula.hard:
            self.add_clause(cl)

    def add_clause(self, clause):
        """
            Add a clause to the solver.
        """

        if clause:
            constr = []
            for l in clause:
                if not abs(l) in self.vmap:
                    v = self.oracle.addVar(vtype=GRB.BINARY, name='x{0}'.format(abs(l)))
                    self.vmap[abs(l)] = v

                constr.append(self.vmap[l] if l > 0 else 1 - self.vmap[-l])

            self.oracle.addConstr(functools.reduce(lambda x, y: x + y, constr) >= 1)
        else:
            # empty clause
            self.oracle.addConstr(0 > 1)

    def compute(self):
        """
            Compute a solution.
        """

        self.oracle.optimize()

        if self.oracle.status != GRB.INFEASIBLE:
            # solution cost
            self.cost = self.oracle.getObjective().getValue()

            # extracting the corresponding model
            self.model = []
            for old, new in self.vmap.items():
                if int(new.X) > 0:
                    self.model.append(old)
                else:
                    self.model.append(-old)

            self.model.sort(key=lambda x: abs(x))
            return self.model

    def enumerate(self):
        """
            Enumerate solutions.
        """

        done = False
        while not done:
            model = self.compute()

            if model != None:
                self.add_clause([-l for l in model])
                yield model
            else:
                done = True

    def delete(self):
        """
            Delete the oracle.
        """

        self.oracle.remove(self.oracle.getConstrs())
        self.oracle.remove(self.oracle.getVars())

        self.oracle = None
