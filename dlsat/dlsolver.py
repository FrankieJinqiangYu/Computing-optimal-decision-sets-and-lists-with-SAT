#!/usr/bin/env python
##
## dlsolve.py
##
##  Created on: Jun 18, 2020
##      Author: Jinqiang Yu
##      E-mail: jinqiang.yu@monash.edu
##

#
#==============================================================================
from pysat.solvers import Glucose3
from pysat.formula import CNF, WCNF
from pysat.examples.rc2 import RC2Stratified


#
#==============================================================================
class DLSolver():
    """
        Class implementing the MaxSAT approach to generating decision lists
    """
    def __init__(self, K, N, data, class_ids=None, LAM=1, mis_weight=1, maxsat=True, sep=False, ml='list', sparse=False):
        self.sep = sep # indicates if a separated model
        self.sparse = sparse # indicates if a sparse model
        self.maxsat = maxsat # indicates if a maxsat model
        self.ml = ml # indicates 'list' or 'hybrid'
        self.K = K  # The number of features.
        self.N = N  # The number of nodes

        # For separated models
        if self.sep:
            self.selected_features = data[:][class_ids] # the feature of data being classified
            self.non_selected_features = [] # the feature of data not being classified


            # For list
            if ml == 'list':
                for ids in range(class_ids + 1, len(data)):
                    self.non_selected_features.extend(data[:][ids]) # data not being classified, not including data that have been classified
            else:
                # for hybrid
                assert ml == 'hybrid'
                for ids in range(0, len(data)):
                    if ids != class_ids:
                        self.non_selected_features.extend(data[:][ids]) # data not being classified, including data that have been classified

            self.M = len(self.selected_features) # The number of items being classified
            self.non_M = len(self.non_selected_features) # The number of item not being classified
            self.C = 1

            # For sparse model, for filtering out data that have been classified
            if sparse:
                num_data = [len(data[:][ids]) for ids in range(class_ids, len(data))]
                self.num_data = [0 for ids in range(class_ids, len(data))]
                for c_id in range(len(num_data)):
                    self.num_data[c_id] = sum(num_data[:c_id + 1])

                self.classified_data = data[:class_ids]
                self.filtered_data = data[class_ids:]

        else:
            # For complete models
            self.selected_features = data[:][0] # the feature of data being classified
            self.M = len(self.selected_features) # The number of items being classified
            self.classes = data[:][1] # the class of data being classified
            self.non_M = 0 # The number of items not being classified
            self.C = 1 if len(self.classes[0]) <= 2 else len(self.classes[0]) #not use one-hot-encoding for binary classes

        self.pi = [[None for j in range(self.K + self.C + 1)] for i in range(self.M + self.non_M + 1)] #pi[i][r]: feature r is true in item i
        self.var_s = [[None for j in range(self.K + self.C + 1)] for i in range(self.N + 1)] # s[j][r]: node j discriminates feature r, var_S[j,K+1] is leaf
        self.var_v = [[None for j in range(self.N + 1)] for i in range(self.M + self.non_M + 1)] # var_v[i][j] is one when item i is valid in node j
        self.var_t = [None for i in range(self.N + 1)] # var_t[j] The value on node j
        self.var_n = [[None for j in range(self.N + 1)] for i in range(self.M + self.non_M + 1)]# n[i][j]  not previously classified by any nodes for the current decision list beforej
        self.var_a0 = [[None for j in range(self.K + self.C + 1)] for i in range(self.N + 1)]# a0[j][r]: the value on node j is 0 and node j discriminates feature r
        self.var_a1 = [[None for j in range(self.K + self.C + 1)] for i in range(self.N + 1)]# a1[j][r]: the value on node j is 1 and node j discriminates feature r
        self.var_aux_vl = [[None for j in range(self.N + 1)] for i in range(self.M + self.non_M + 1)] # vl[i][j]: node j is a leaf and item i is valid at node j
        self.var_u = [None for i in range(self.N + 1)]  # u[j] node j is unused

        #if self.ml == 'hybrid'
        self.var_x = [None for i in range(self.N + 1)] # x[j]: node j ends a decision lis

        # For sparse models
        if sparse:
            self.var_m = [None for i in range(self.M + self.non_M + 1)] # m[i] item i is misclassified
            self.LAM = LAM
            self.mis_weight = mis_weight

        # cnf
        if self.maxsat:
            self.cnf = WCNF()
        else:
            self.cnf = CNF()

        self.solver = None
        self.runtime = 0

        self.var2ids = {} #{variable: index in sat}
        self.ids2var = {} #{index in sat: variable}
        self.cnf_ids = 0 # record the variable index in cnf
        self.all_var_sol = None # store a solution
        self.basic_keys = ["s", "t", "x", "u", "m"] # basic keys for building decision lists

        #Parse values to pi
        if self.sep:
            for i in range(1, self.M + self.non_M + 1):
                for r in range(1, self.K + 1):
                    if i < (self.M + 1):
                        self.pi[i][r] = self.selected_features[i - 1][r - 1]
                    else:
                        self.pi[i][r] = self.non_selected_features[i - self.M - 1][r - 1]
        else:
            for i in range(1, self.M + self.non_M + 1):
                for r in range(1, self.K + 1):
                    self.pi[i][r] = self.selected_features[i - 1][r - 1]

            for i in range(len(self.classes)):
                for j in range(self.C):
                    if self.classes[i][j] == 1:
                        self.pi[i + 1][self.K + j + 1] = 1
                    else:
                        self.pi[i + 1][self.K + j + 1] = 0

        # Assign the name of a variable
        def add_vars_1d(vars, var_name):
            for i in range(0, len(vars)):
                vars[i] = var_name + "[" + str(i) + "]"

        def add_vars_2d(vars, var_name):
            for i in range(0, len(vars)):
                for j in range(0, len(vars[0])):
                    vars[i][j] = var_name + "[" + str(i) + "," + str(j) + "]"

        add_vars_2d(self.var_v, "v")
        add_vars_2d(self.var_s, "s")
        add_vars_1d(self.var_t, "t")
        add_vars_2d(self.var_a0, "a0")
        add_vars_2d(self.var_a1, "a1")
        add_vars_2d(self.var_aux_vl, "aux_vl")
        add_vars_2d(self.var_n, "n")
        add_vars_1d(self.var_x, "x")
        add_vars_1d(self.var_u, "u")
        if sparse:
            add_vars_1d(self.var_m, "m")


    def get_cnf_id(self, var, value=1):
        """

            Input a variable, then get the index of the variable in cnf
        """
        if var not in self.var2ids:
            self.cnf_ids += 1
            self.var2ids[var] = self.cnf_ids

            assert self.cnf_ids not in self.ids2var
            self.ids2var[self.cnf_ids] = var

        return self.var2ids[var] if value == 1 else -self.var2ids[var]


    def encode_constraints(self):
        """

            Encode the constraints
        """

        # A node uses only one feature (or the class feature)
        # Exactly one u[j] + S[j1] + S[j2]... + S[jf] + S[jc1] +  S[jc2] + ....
        def cons_1(self):
            for j in range(1, self.N + 1):
                or_f = []
                or_f.append(self.get_cnf_id(self.var_u[j]))
                for f in range(1, self.K + self.C + 1):
                    or_f.append(self.get_cnf_id(self.var_s[j][f]))
                self.cnf.append(or_f)

                for a in range(len(or_f)):
                    for b in range(a + 1, len(or_f)):
                        new_f = []
                        new_f.append(-or_f[a])
                        new_f.append(-or_f[b])
                        self.cnf.append(new_f)

        #
        # Auxiliary: a0[j, r] <-> s[j, r] /\ -t[j], a1[j, r] <-> s[j, r] /\ t[j]
        def cons_au(self):
            for j in range(1, self.N + 1):
                for r in range(1, self.K + 1):
                    # Auxiliary: a0[j, r] <-> s[j, r] /\ -t[j]
                    f1 = []
                    f1.append(self.get_cnf_id(self.var_a0[j][r], '-'))
                    f1.append(self.get_cnf_id(self.var_s[j][r]))
                    self.cnf.append(f1)

                    f2 = []
                    f2.append(self.get_cnf_id(self.var_a0[j][r], '-'))
                    f2.append(self.get_cnf_id(self.var_t[j], '-'))
                    self.cnf.append(f2)

                    f3 = []
                    f3.append(self.get_cnf_id(self.var_a0[j][r]))
                    f3.append(self.get_cnf_id(self.var_s[j][r], '-'))
                    f3.append(self.get_cnf_id(self.var_t[j]))
                    self.cnf.append(f3)

                    # Auxiliary: a1[j, r] <-> s[j, r] /\ t[j]
                    f1 = []
                    f1.append(self.get_cnf_id(self.var_a1[j][r], '-'))
                    f1.append(self.get_cnf_id(self.var_s[j][r]))
                    self.cnf.append(f1)

                    f2 = []
                    f2.append(self.get_cnf_id(self.var_a1[j][r], '-'))
                    f2.append(self.get_cnf_id(self.var_t[j]))
                    self.cnf.append(f2)

                    f3 = []
                    f3.append(self.get_cnf_id(self.var_a1[j][r]))
                    f3.append(self.get_cnf_id(self.var_s[j][r], '-'))
                    f3.append(self.get_cnf_id(self.var_t[j], '-'))
                    self.cnf.append(f3)

        # If a nodejis unused then so are all the following nodes
        # u[j] -> u[j+1]
        def cons_2(self):
            for j in range(1, self.N):
                f1 = []
                f1.append(self.get_cnf_id(self.var_u[j], '-'))
                f1.append(self.get_cnf_id(self.var_u[j + 1]))
                self.cnf.append(f1)

        # The last used node is a leaf (cons_3, 4)
        # u[j+1] -> u[j]  \/  {for any c in C} s[j,c]
        def cons_3(self):
            for j in range(1, self.N):
                f1 = []
                f1.append(self.get_cnf_id(self.var_u[j+1], '-'))
                f1.append(self.get_cnf_id(self.var_u[j]))
                for c in range(self.K + 1, self.K + self.C + 1):
                    f1.append(self.get_cnf_id(self.var_s[j][c]))
                self.cnf.append(f1)

        # u[j+1] -> u[N]  \/  (for any c in C, s[j,c])
        def cons_4(self):
            for j in range(1, self.N):
                f1 = []
                f1.append(self.get_cnf_id(self.var_u[j + 1], '-'))
                f1.append(self.get_cnf_id(self.var_u[self.N]))
                for c in range(self.K + 1, self.K + self.C + 1):
                    f1.append(self.get_cnf_id(self.var_s[self.N][c]))
                self.cnf.append(f1)

        # all examples are not previously classified at the first node
        # n[i, 1] for all data
        def cons_5(self):
            for i in range(1, self.M + self.non_M + 1):
                self.cnf.append([self.get_cnf_id(self.var_n[i][1])])

        # An example e_i is previously unclassified at node j+ 1
        # iff it was previously unclassified, and either j is a not leaf node or it was invalid at the previous leaf node
        # (so not classified by the rule that finished there)
        # n[i, j+1] <->  n[i,j] /\ ((for and c in C, -s[j,c]) \/ -v[i,j])
        def cons_6(self):
            for i in range(1, self.M + self.non_M + 1):
                for j in range(1, self.N):
                    f1 = []
                    f1.append(self.get_cnf_id(self.var_n[i][j + 1], '-'))
                    f1.append(self.get_cnf_id(self.var_n[i][j]))
                    self.cnf.append(f1)

                    f2 = []
                    f2.append(self.get_cnf_id(self.var_n[i][j + 1], '-'))
                    f2.append(self.get_cnf_id(self.var_v[i][j], '-'))
                    for c in range(self.K + 1, self.K + self.C + 1):
                        or_f = f2[:]
                        or_f.append(self.get_cnf_id(self.var_s[j][c], '-'))
                        self.cnf.append(or_f)

                    f3 = []
                    f3.append(self.get_cnf_id(self.var_n[i][j + 1]))
                    f3.append(self.get_cnf_id(self.var_n[i][j], '-'))
                    for c in range(self.K + 1, self.K + self.C + 1):
                        f3.append(self.get_cnf_id(self.var_s[j][c]))
                    self.cnf.append(f3)

                    f4 = []
                    f4.append(self.get_cnf_id(self.var_n[i][j + 1]))
                    f4.append(self.get_cnf_id(self.var_n[i][j], '-'))
                    f4.append(self.get_cnf_id(self.var_v[i][j]))
                    self.cnf.append(f4)

        # All examples are valid at the first node
        # v[i,1]
        def cons_7(self):
            for i in range(1, self.M + self.non_M + 1):
                cnf_id = [self.get_cnf_id(self.var_v[i][1])]
                self.cnf.append(cnf_id)

        # An example e_i is valid at nodej+ 1 iff j is a leaf node and it was previously unclassified,
        # or e_i is valid at node j and e_i and node j agree on the value of the feature s_jr selected for that node
        # v[i, j+1] <-> (s[j,c] /\ n[i,j+1]) \/ (v[i,j] /\ {for any r in K} (s[j,r] /\ (t[j] = pi_[i, r])))
        def cons_8(self):
            for i in range(1, self.M + self.non_M + 1):
                for j in range(1, self.N):
                    f1 = []
                    f1.append(self.get_cnf_id(self.var_v[i][j + 1], '-'))
                    f1.append(self.get_cnf_id(self.var_v[i][j]))
                    for c in range(self.K+1, self.K+self.C+1):
                        f1.append(self.get_cnf_id(self.var_s[j][c]))
                    self.cnf.append(f1)

                    f2 = []
                    f2.append(self.get_cnf_id(self.var_v[i][j + 1], '-'))
                    f2.append(self.get_cnf_id(self.var_v[i][j]))
                    f2.append(self.get_cnf_id(self.var_n[i][j + 1]))
                    self.cnf.append(f2)

                    f3 = []
                    f3.append(self.get_cnf_id(self.var_v[i][j + 1], '-'))
                    for r in range(1, self.K + 1):
                        if self.pi[i][r] == 0:
                            f3.append(self.get_cnf_id(self.var_a0[j][r]))
                        else:
                            assert self.pi[i][r] == 1
                            f3.append(self.get_cnf_id(self.var_a1[j][r]))
                    for c in range(self.K+1, self.K+self.C+1):
                        f3.append(self.get_cnf_id(self.var_s[j][c]))
                    self.cnf.append(f3)

                    f4 = []
                    f4.append(self.get_cnf_id(self.var_v[i][j + 1], '-'))
                    for r in range(1, self.K + 1):
                        if self.pi[i][r] == 0:
                            f4.append(self.get_cnf_id(self.var_a0[j][r]))
                        else:
                            assert self.pi[i][r] == 1
                            f4.append(self.get_cnf_id(self.var_a1[j][r]))
                    f4.append(self.get_cnf_id(self.var_n[i][j + 1]))
                    self.cnf.append(f4)

                    or_f = []
                    or_f.append(self.get_cnf_id(self.var_v[i][j + 1]))
                    or_f.append(self.get_cnf_id(self.var_n[i][j + 1], '-'))
                    for c in range(self.K + 1, self.K + self.C + 1):
                        f5 = or_f[:]
                        f5.append(self.get_cnf_id(self.var_s[j][c], '-'))
                        self.cnf.append(f5)

                    f6 = []
                    f6.append(self.get_cnf_id(self.var_v[i][j + 1]))
                    f6.append(self.get_cnf_id(self.var_v[i][j], '-'))
                    for r in range(1, self.K + 1):
                        ff = f6[:]
                        if self.pi[i][r] == 0:
                            ff.append(self.get_cnf_id(self.var_a0[j][r], '-'))
                        else:
                            assert self.pi[i][r] == 1
                            ff.append(self.get_cnf_id(self.var_a1[j][r], '-'))
                        self.cnf.append(ff)

        # If example e_i is valid at a leaf node j, they should agree on the class feature
        # s[j,c] /\ v[i, j] -> (t[j] = c[i])
        def cons_9(self):
            if self.C == 1:
                for i in range(1, self.M + 1):
                    for j in range(1, self.N + 1):
                        f1 = []
                        f1.append(self.get_cnf_id(self.var_s[j][self.K + 1], '-'))
                        f1.append(self.get_cnf_id(self.var_v[i][j], '-'))
                        if self.pi[i][self.K + 1] == 1:
                            f1.append(self.get_cnf_id(self.var_t[j]))
                        else:
                            f1.append(self.get_cnf_id(self.var_t[j], '-'))
                        self.cnf.append(f1)
            else:
                assert self.C >= 3
                # s[j,c] /\ v[i, j] -> c[i]
                for i in range(1, self.M + 1):
                    for j in range(1, self.N + 1):
                        for c in range(self.K+1, self.K+self.C+1):
                            if self.pi[i][c] != 1:
                                f1 = []
                                f1.append(self.get_cnf_id(self.var_s[j][c], '-'))
                                f1.append(self.get_cnf_id(self.var_v[i][j], '-'))
                                self.cnf.append(f1)

        # When there are 3 or more classes we restrict leaf nodes to only consider
        # true examples of the class
        # s[j,c] -> t[j]
        def cons_10(self):
            if self.C >= 3:
                for j in range(1, self.N + 1):
                    for c in range(self.K+1, self.K+self.C+1):
                        f1 = []
                        f1.append(self.get_cnf_id(self.var_s[j][c], '-'))
                        f1.append(self.get_cnf_id(self.var_t[j]))
                        self.cnf.append(f1)

        # For every example there should be at least one leaf node where it is valid:
        # {For any j in N} ( {for any c in C} s[j,c] /\ v[i,j] )
        def cons_11(self):
            # Auxiliary variables vl[i,j] <-> {for any c in C }s[jc] /\ v[ij]
            for i in range(1, self.M + self.non_M + 1):
                for j in range(1, self.N + 1):
                    f1 = []
                    f1.append(self.get_cnf_id(self.var_aux_vl[i][j], '-'))
                    for c in range(self.K+1, self.K+self.C+1):
                        f1.append(self.get_cnf_id(self.var_s[j][c]))
                    self.cnf.append(f1)

                    f2 = []
                    f2.append(self.get_cnf_id(self.var_aux_vl[i][j], '-'))
                    f2.append(self.get_cnf_id(self.var_v[i][j]))
                    self.cnf.append(f2)

                    or_f = []
                    or_f.append(self.get_cnf_id(self.var_aux_vl[i][j]))
                    or_f.append(self.get_cnf_id(self.var_v[i][j], '-'))
                    for c in range(self.K+1, self.K+self.C+1):
                        f3 = or_f[:]
                        f3.append(self.get_cnf_id(self.var_s[j][c], '-'))
                        self.cnf.append(f3)

            # For selected items, at least valid at one leaf
            # {For all i in selected items} {For any j in N} (vl[i, j])
            for i in range(1, self.M + 1):
                or_f = []
                for j in range(1, self.N + 1):
                    or_f.append(self.get_cnf_id(self.var_aux_vl[i][j]))
                self.cnf.append(or_f)

            # For selected items, not valid at all leaves
            # For non-selected items, {For any j in N} ( {for any c in C} s[j,c] /\ v[i,j] )
            # - ({for any c in C} s[jc] /\ v[i, j])
            # (-vl[i,j])
            if self.sep:
                for i in range(self.M + 1, self.M + self.non_M + 1):
                    for j in range(1, self.N + 1):
                        self.cnf.append([self.get_cnf_id(self.var_aux_vl[i][j], '-')])


        # The las used node ends a decision list
        # u[j+1] /\ -u[j] -> x[j]
        def cons_12(self):
            for j in range(1, self.N):
                f1 = []
                f1.append(self.get_cnf_id(self.var_u[j+1], '-'))
                f1.append(self.get_cnf_id(self.var_u[j]))
                f1.append(self.get_cnf_id(self.var_x[j]))
                self.cnf.append(f1)

        # An end node is always a leaf
        # x[j] -> s[j,c1] \/ s[j,c2] ....
        def cons_13(self):
            for j in range(1, self.N + 1):
                f1 = []
                f1.append(self.get_cnf_id(self.var_x[j], '-'))
                for c in range(self.K+1, self.K+self.C+1):
                    f1.append(self.get_cnf_id(self.var_s[j, c]))
                self.cnf.append(f1)

        # An  example e_i is previously unclassified at node j+ 1
        # iff j is an end of decision list node, or it was previously unclassified, and either j is a not leaf node
        # or it was invalid at the previous leaf node (so not classified by the rule that finished there)
        # n[i,j+1] <->x[j] \/ (n[i,j] /\( ({any c in C}-s[j,c]) \/ -v[i,j] ) )
        def cons_14(self):
            for i in range(1, self.M + self.non_M + 1):
                for j in range(1, self.N):
                    f1 = []
                    f1.append(self.get_cnf_id(self.var_x[j], '-'))
                    f1.append(self.get_cnf_id(self.var_n[i][j + 1]))
                    self.cnf.append(f1)

                    f2 = []
                    f2.append(self.get_cnf_id(self.var_n[i][j], '-'))
                    f2.append(self.get_cnf_id(self.var_n[i][j + 1]))
                    for c in range(self.K+1, self.K+self.C+1):
                        f2.append(self.get_cnf_id(self.var_s[j, c]))
                    self.cnf.append(f2)

                    f3 = []
                    f3.append(self.get_cnf_id(self.var_v[i][j]))
                    f3.append(self.get_cnf_id(self.var_n[i][j], '-'))
                    f3.append(self.get_cnf_id(self.var_n[i][j + 1]))
                    self.cnf.append(f3)

                    f4 = []
                    f4.append(self.get_cnf_id(self.var_n[i][j]))
                    f4.append(self.get_cnf_id(self.var_x[j]))
                    f4.append(self.get_cnf_id(self.var_n[i][j + 1], '-'))
                    self.cnf.append(f4)

                    f5 = []
                    f5.append(self.get_cnf_id(self.var_v[i][j], '-'))
                    f5.append(self.get_cnf_id(self.var_x[j]))
                    f5.append(self.get_cnf_id(self.var_n[i][j + 1], '-'))
                    for c in range(self.K+1, self.K+self.C+1):
                        or_f = f5[:]
                        or_f.append(self.get_cnf_id(self.var_s[j, c], '-'))
                        self.cnf.append(or_f)

        # If example e_i is valid at a leaf node j
        # then they agree on the class feature or the item is misclassified
        # s[j,c] /\ v[i, j] -> t[j] = c[i] \/ m[i]
        def cons_15(self):
            if self.C == 1:
                for i in range(1, self.M + 1):
                    for j in range(1, self.N + 1):
                        f1 = []
                        f1.append(self.get_cnf_id(self.var_s[j][self.K + 1], '-'))
                        f1.append(self.get_cnf_id(self.var_v[i][j], '-'))
                        if self.pi[i][self.K + 1] == 1:
                            f1.append(self.get_cnf_id(self.var_t[j]))
                        else:
                            f1.append(self.get_cnf_id(self.var_t[j], '-'))
                        f1.append(self.get_cnf_id(self.var_m[i]))
                        self.cnf.append(f1)
            else:
                assert self.C >= 3
                # s[j,c] /\ v[i, j] -> c[i] \/ m[i]
                for i in range(1, self.M + 1):
                    for j in range(1, self.N + 1):
                        for c in range(self.K + 1, self.K + self.C + 1):
                            if self.pi[i][c] != 1:
                                f1 = []
                                f1.append(self.get_cnf_id(self.var_s[j][c], '-'))
                                f1.append(self.get_cnf_id(self.var_v[i][j], '-'))
                                f1.append(self.get_cnf_id(self.var_m[i]))
                                self.cnf.append(f1)

        # For every example there should be at least one leaf node where it is valid:
        # m[i] \/ {For any j in N} ( {for any c in C} s[j,c] /\ v[i,j] )
        def cons_16(self):
            # Auxiliary variables vl[i,j] <-> {for any c in C }s[jc] /\ v[ij]
            for i in range(1, self.M + self.non_M + 1):
                for j in range(1, self.N + 1):
                    f1 = []
                    f1.append(self.get_cnf_id(self.var_aux_vl[i][j], '-'))
                    for c in range(self.K + 1, self.K + self.C + 1):
                        f1.append(self.get_cnf_id(self.var_s[j][c]))
                    self.cnf.append(f1)

                    f2 = []
                    f2.append(self.get_cnf_id(self.var_aux_vl[i][j], '-'))
                    f2.append(self.get_cnf_id(self.var_v[i][j]))
                    self.cnf.append(f2)

                    or_f = []
                    or_f.append(self.get_cnf_id(self.var_aux_vl[i][j]))
                    or_f.append(self.get_cnf_id(self.var_v[i][j], '-'))
                    for c in range(self.K + 1, self.K + self.C + 1):
                        f3 = or_f[:]
                        f3.append(self.get_cnf_id(self.var_s[j][c], '-'))
                        self.cnf.append(f3)

            # For selected items, at least valid at one leaf or misclassified
            # {For all i in selected items} m[i] \/ {For any j in N} (vl[i, j])
            for i in range(1, self.M + 1):
                or_f = []
                or_f.append(self.get_cnf_id(self.var_m[i]))
                for j in range(1, self.N + 1):
                    or_f.append(self.get_cnf_id(self.var_aux_vl[i][j]))
                self.cnf.append(or_f)

            # For selected items, not valid at all leaves or misclassified
            # For non-selected items, {For any j in N} ( {for any c in C} s[j,c] /\ v[i,j] )
            # mi \/ - ({for any c in C} s[jc] /\ v[i, j])
            # mi \/ (-vl[i,j])
            if self.sep:
                for i in range(self.M + 1, self.M + self.non_M + 1):
                    for j in range(1, self.N + 1):
                        f1 = []
                        f1.append(self.get_cnf_id(self.var_m[i]))
                        f1.append(self.get_cnf_id(self.var_aux_vl[i][j], '-'))
                        self.cnf.append(f1)

        #####################

        # if it is a hybrid model
        if self.ml == 'hybrid' and not self.sep:
                cons_1(self)
                cons_au(self)
                cons_2(self)
                cons_3(self)
                cons_4(self)
                cons_5(self)
                cons_7(self)
                cons_8(self)
                cons_9(self)
                cons_10(self)
                cons_11(self)
                cons_12(self)
                cons_13(self)
                cons_14(self)
                # Maximise u[j]
                for j in range(1, self.N + 1):
                    self.cnf.append([self.get_cnf_id(self.var_u[j])], weight=1)

        else:
            # if it is a perfect model
            if self.maxsat and not self.sparse:
                cons_1(self)
                cons_au(self)
                cons_2(self)
                cons_3(self)
                cons_4(self)
                cons_5(self)
                cons_6(self)
                cons_7(self)
                cons_8(self)
                if not self.sep:
                    cons_9(self)
                cons_10(self)
                cons_11(self)
                # Maximise u[j]
                for j in range(1, self.N + 1):
                    self.cnf.append([self.get_cnf_id(self.var_u[j])], weight=1)

            # if it is a sparse model
            if self.maxsat and self.sparse:
                cons_1(self)
                cons_au(self)
                cons_2(self)
                cons_3(self)
                cons_4(self)
                cons_5(self)
                cons_6(self)
                cons_7(self)
                cons_8(self)
                cons_10(self)
                if not self.sep:
                    cons_15(self)
                cons_16(self)

                # objective function
                for i in range(1, self.M + self.non_M + 1):
                    self.cnf.append([self.get_cnf_id(self.var_m[i], '-')], weight=self.mis_weight)

                for j in range(1, self.N + 1):
                    self.cnf.append([self.get_cnf_id(self.var_u[j])], weight=self.LAM)




    def solve(self):
        """

            Compute an optimal decision list
        """
        self.encode_constraints()
        self.all_var_sol = []
        ds_var_sol = []

        if self.maxsat:
            self.solver = RC2Stratified(self.cnf, solver='g3', adapt=True, exhaust=True, incr=False, minz=True, trim=0, verbose=0)
            s = self.solver.compute()
            self.runtime = self.solver.oracle_time()
            if s is None:
                return None
            else:
                # Sparse variables, e.g ds_var_sol: [[(s[1,1], True), (s[1,2], False)]]
                for ids in s:
                    var_name = self.ids2var[abs(ids)]
                    var_val = 1 if ids > 0 else 0
                    self.all_var_sol.append((var_name, var_val))
                    if var_name[0] in self.basic_keys:
                        ds_var_sol.append((var_name, var_val))

                return [ds_var_sol]

        else:
            self.solver = Glucose3(bootstrap_with=self.cnf, use_timer=True)
            found = self.solver.solve()
            self.runtime = self.solver.time()
            if not found:
                return None
            else:
                s = self.solver.get_model()
                # Sparse variables, e.g ds_var_sol: [[(s[1,1], True), (s[1,2], False)]]
                for ids in s:
                    var_name = self.ids2var[abs(ids)]
                    var_val = 1 if ids > 0 else 0
                    self.all_var_sol.append((var_name, var_val))
                    if var_name[0] in self.basic_keys:
                        ds_var_sol.append((var_name, var_val))

                return [ds_var_sol]






