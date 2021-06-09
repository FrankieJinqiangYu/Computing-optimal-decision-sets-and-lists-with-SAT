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


#
#==============================================================================
class DL():
    """

        Class of decisoin lists
    """
    def __init__(self, K, N, C=1):
        self.K = K #The number of features.
        self.N = N #The number of nodes
        self.C = 1 if C==2 else C # The number of classes. If it is a binary class, c=1. If classes > 3, using one-hot-encoding
        self.valid_N = 0 # The number of used nodes

        self.s = [[False for j in range(self.K + self.C + 1)] for i in range(self.N + 1)] # S[i][f]: Node i discriminates on feature f
        self.t = [False for i in range(self.N + 1)] # t[j] The value on node j
        self.u = [False for i in range(self.N + 1)] # u[j] Node j is unused
        self.x = [False for i in range(self.N + 1)] # node j ends a decision list

        self.mis_count = 0 # Misclassification number
        self.node = [] # Solution. e.g [('f1', True),('class1', True)]


    def get_num_literals(self):
        """

            Get the number of literals in decision lists, not including rules
        """
        num_rule = 0
        for j in range(1, self.N+1):
            for c in range(self.K+1, self.K+self.C+1):
                if self.s[j][c]:
                    num_rule += 1
                    break
        return self.valid_N  - num_rule


    def parse_solution(self, sols, num_nodes=None):
        """

            Parse solution to generate list
            solution = [(key, value),(key, value)....] key = (''v[1], 'v[2]', ....'v[N+1]',c'[1]'....) and other variables keys
            num_nodes: A list contains the number of used nodes in each class, e.g. [2,3] 2 nodes for class
        """
        if num_nodes is not None:
            self.valid_N = sum(num_nodes)
        var_name_lookup = {"s": self.s, "t": self.t, "u": self.u, "m": None, "x": self.x}
        sol_ids = 0
        num_pre_nodes = 0
        for sol in sols:
            for var, value in sol:
                var_name = var[: var.find("[")].strip()
                var_index = var[var.find("[") + 1: var.find("]")].strip().split(",")
                assert var_name in var_name_lookup
                if var_name == "t" or var_name == "u" or var_name == "x":
                    var_name_lookup[var_name][int(var_index[0]) + num_pre_nodes] = True if value == 1 else False
                elif var_name != "m":
                    assert var_name == "s"
                    j = int(var_index[0])
                    f = int(var_index[1])
                    var_name_lookup[var_name][j + num_pre_nodes][f] = True if value == 1 else False

            if num_nodes is None:
                for j in range(1, self.N + 1):
                    if not self.u[j]:
                      self.valid_N += 1
            else:
                num_pre_nodes += num_nodes[sol_ids]
                sol_ids += 1



    def generate_node(self, num_nodes=None, cur_class=None, class_order=None):
        """

            Generate nodes using variables
            Example: [('f1', True), ('f2', False), ('class1', True)]
        """
        for j in range(1, self.N + 1):
            if not self.u[j]:
                value = self.t[j]

                if sum(self.s[j][self.K+1: ]) == 1:
                    feature = "class"
                    for c in range(self.K+1, self.K+self.C+1):
                        if self.s[j][c]:
                            c_ids = c - self.K
                            break
                else:
                    feature = None
                    for f in range(1, self.K + 1):
                        if self.s[j][f]:
                            feature = "f" + str(f)
                            break

                if feature is not None:
                    if feature != 'class':
                        self.node.append((feature, value))
                    else:
                        if num_nodes is None:
                            if cur_class is None:
                                if self.C >= 3:
                                    self.node.append((feature + str(c_ids), value))
                                else:
                                    self.node.append((feature + '1' if value else feature + '2', value))
                            else:
                                self.node.append((feature+str(cur_class), True))
                        else:
                            num = 0
                            for c_ids in range(len(num_nodes)):
                                num += num_nodes[c_ids]
                                if len(self.node) < num:
                                    self.node.append((feature + str(class_order[c_ids]), True))
                                    break


    def validate(self, data_features, data_classes, num_nodes, num_class=2, ml='list'):

        """

            Validation, calculate the accuracy. for separated models
            num_nodes: A list contains the number of used and unused nodes in each class
            num_class: the number of classes in datasets
            c_order: the class order in the separated decision lists
        """
        self.mis_count = 0 # the number of miclassified items
        self.lit_used = [] # the literals used

        # Assign a default rule if there was no rules in a decision list.
        # The rule is the class with the most number of training datasets.
        if self.node == []:
            class_rec = dict() #Record the number if item in each class
            for c in data_classes:
                for i in range(num_class):
                    if c[i]:
                        class_rec[i+1] = class_rec.get(i+1, 0) + 1
                        break

            max_class_ids = 1 # The class index with most items
            max_class = 0 # The most number of items in the max class

            for i in range(1, num_class+1):
                if class_rec.get(i, 0) > max_class:
                    max_class_ids = i
                    max_class = class_rec.get(i, 0)
            self.node.append(('class' + str(max_class_ids), True))
            self.valid_N = 1

        # The calculations are different for lists and hybrid models (set of lists)
        # For lists
        if ml == 'list':
            for f, c in zip(data_features, data_classes):
                pre_match = True # Item is valid at the current node
                pre_leaf = False # The previous node is a leaf (rule)
                mis = False # True if the item is misclassified
                found = False # Item is found
                lit_used = 1

                # Go through all nodes in the decision list
                for node in self.node:
                    node_value = node[1] # The value of a feature, True or False
                    if node[0].lower().startswith('f'): # If the current node is a literal
                        node_f = int(node[0][1:]) # The index of the feature, e.g. 1 (the first)
                        if pre_leaf:
                            pre_leaf = False
                            if f[node_f - 1] == node_value: # If item value == node value
                                pre_match = True # For the next node
                            else:
                                pre_match = False
                        elif pre_match and f[node_f - 1] == node_value: # pre_match is unchanged, pre_match = True
                            pass
                        else:
                            # Item is not valid at the current node OR the feature value is different for the item and the list
                            assert f[node_f - 1] != node_value or not pre_match
                            pre_match = False # Item is not valid at the next node
                    else:
                        # If the current node is a leaf (rule)
                        assert node[0][ : 5].lower().strip() == 'class'
                        pre_leaf = True
                        item_class = list(c).index(1)
                        # If the item is valid at the node but the classes are different
                        # The item is misclassified
                        if pre_match:
                            found = True
                            if str(item_class + 1) != node[0][ 5: ]:
                                mis = True
                            break
                        pre_match = True

                    lit_used += 1

                if not found:
                    mis = True
                if mis:
                    self.mis_count += 1
                    self.lit_used.append(len(self.node))
                else:
                    self.lit_used.append(lit_used)

        elif ml == 'hybrid':
            # For hybrid
            if sum(self.x) == 0:
                # Divide nodes based on the class
                self.set_nodes = {c+1: None for c in range(num_class)}

                if self.valid_N != 1:
                    self.set_nodes[1] = self.node[: num_nodes[0]]
                    self.set_nodes[num_class] = self.node[sum(num_nodes[0 : num_class-1]): ]

                    if num_class > 2:
                        for c in range(1, num_class-1):
                            self.set_nodes[c+1] = self.node[sum(num_nodes[0 : c]): sum(num_nodes[0: c+1])]
                else:
                    assert len(self.node) == 1
                    node_class = int(self.node[0][0][5:])
                    for c in range(1, num_class+1):
                        if c == node_class:
                            self.set_nodes[c] = self.node[:]
                        else:
                            self.set_nodes[c] = []

                data_predicted = [0 for num_data in range(len(data_features))] # Record classification results

                # Go through all items
                for data_ids in range(len(data_features)):
                    lit_used_list = {}
                    f = data_features[data_ids] # Item features
                    c = data_classes[data_ids] # Item class
                    for iii in range(len(c)):
                        if c[iii]:
                            data_class = iii+1
                            break
                    mis = False

                    #Go through all lists in the set
                    for c_ids in range(1, num_class+1):
                        lit_used = 0
                        if mis:
                            break
                        found = False # Item is found
                        pre_match = True # Item is valid at the current node
                        pre_leaf = False # Previous node is a leaf node
                        nodes = self.set_nodes[c_ids][:] # The decision list for the class

                        # If the decision list is empty
                        if nodes == []:

                            # If item class is as same as the class in the decision list
                            # The item is misclassified
                            if c[c_ids-1] == 1:
                                mis = True
                                break
                        else:
                            #Go through all literals in the decision list
                            for node in nodes:
                                node_value = node[1] # The value for the feature, True or False
                                lit_used += 1

                                if node[0].lower().startswith('f'): # If the node is a literal not a leaf node
                                    node_f = int(node[0][1:]) # The index of the feature
                                    if pre_leaf:
                                        pre_leaf = False
                                        if f[node_f - 1] == node_value: # The feature values are the same for the item and the decisoin list
                                            pre_match = True # For the next node
                                        else:
                                            pre_match = False
                                    elif pre_match and f[node_f - 1] == node_value: # pre_match is unchanged, pre_match = True
                                        pass
                                    else:
                                        pre_match = False
                                else:
                                    # If the node is a leaf node (rule)
                                    assert node[0][ : 5].lower().strip() == 'class'
                                    pre_leaf = True
                                    item_class = list(c).index(1)
                                    if pre_match:
                                        found = True
                                        lit_used_list[data_class] = lit_used

                                        # If the item is valid at the node but the classes are different
                                        # The item is misclassified
                                        if str(item_class + 1) != node[0][ 5: ]:
                                            mis = True
                                        break

                                    pre_match = True

                            if found and mis:
                                break

                            # The classes in the item and decision list are the same, but the item is not classified
                            # The item is misclassified
                            if not found and str(item_class + 1) == node[0][ 5: ]:
                                mis = True
                                break
                    # Record misclassification
                    if mis:
                        data_predicted[data_ids] = 1
                        self.lit_used.append(len(self.node))
                    else:
                        self.lit_used.append(lit_used_list[data_class])

                # Count the misclssified items

                self.mis_count = sum(data_predicted)


    def get_miss(self, num_items, sol):
        """

            Get the misclassified items
            Used for filtering out classified and misclassified items
        """
        m = [0 for i in range(num_items + 1)]
        var_name_lookup = {"m": m}
        misclas_ids = []
        for var, value in sol:
            var_name = var[: var.find("[")].strip()
            var_index = var[var.find("[") + 1: var.find("]")].strip().split(",")
            if var_name == "m":
                var_name_lookup[var_name][int(var_index[0])] = value
                if value == 1:
                    misclas_ids.append(int(var_index[0]))

        # Return variable m and the indices of misclassified items
        return sum(m[1:]), misclas_ids


    def generate_list(self, feature_dict, feature_names, ml='list'):
        """

            Generate a decision list
        """
        self.tuple_list = [() for i in range(len(self.node))]
        for node_ids in range(len(self.node)):
            node_feature, node_value = self.node[node_ids]
            node_value = 1 if node_value is True else 0
            is_class = True if node_feature.strip().lower().startswith('class') else False
            if not is_class:
                assert node_feature.strip().lower().startswith('f'), 'not a feature'
                node_feature_ids = int(node_feature[1:]) - 1
            else:
                node_feature_ids = int(node_feature[5:]) - 1

            self.tuple_list[node_ids] = (feature_names[node_feature_ids] if not is_class else feature_names[-1],
                                         feature_dict[feature_names[node_feature_ids]][node_value] if not is_class
                                         else feature_dict[feature_names[-1]][node_feature_ids])
        # Node to list
        if ml == 'list':
            if len(self.node) == 0:
                self.list = None
            else:
                if len(self.node) == 1 or self.node[0][0].lower().startswith('class'):

                    self.list = self.tuple_list[0][0] + '=' + str(self.tuple_list[0][1])
                else:
                    self.list = 'IF ' + self.tuple_list[0][0] + '=' + str(self.tuple_list[0][1])
                    for ins in range(1, len(self.node)):
                        pre = self.node[ins - 1]
                        curr = self.node[ins]
                        if pre[0].lower().startswith('f') and curr[0].lower().startswith('f'):
                            self.list = self.list + ' AND ' + self.tuple_list[ins][0] + '=' + str(self.tuple_list[ins][1])
                        elif  pre[0].lower().startswith('f') and curr[0].lower().startswith('class'):
                            self.list = self.list + ' THEN ' + self.tuple_list[ins][0] + '=' + str(self.tuple_list[ins][1]) + '\n'
                        elif pre[0].lower().startswith('class') and curr[0].lower().startswith('f'):
                            self.list = self.list + 'ELIF ' + self.tuple_list[ins][0] + '=' + str(self.tuple_list[ins][1])
                        elif pre[0].lower().startswith('class') and curr[0].lower().startswith('class'):
                            self.list = self.list + 'ELSE ' + self.tuple_list[ins][0] + '=' + str(self.tuple_list[ins][1]) + '\n'

