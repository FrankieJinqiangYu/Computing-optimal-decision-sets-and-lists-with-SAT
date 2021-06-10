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
from __future__ import print_function
import sys
import getopt
import os
from dlsolver import DLSolver
from dataprocessing import data_processing
from dl import DL
import math
import itertools
import random


#
#==============================================================================

def parse_options():
    """
        Parses command-line options.
    """

    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   'a:h:l:m:n:o:',
                                   ['help',
                                    'lam==',
                                    'lambda==',
                                    'maxsat',
                                    'mxsat',
                                    'ml',
                                    'model==',
                                    'node==',
                                    'order==',
                                    'sep',
                                    'separate',
                                    'spa',
                                    'sparse'])
    except getopt.GetoptError as err:
        sys.stderr.write(str(err).capitalize() + '\n')
        usage()
        sys.exit(1)

    train_filename = args[0]

    lam = 0.005
    maxsat = True # If a maxsat model
    ml = 'list'  # 'list, 'hybrid'
    N = 1
    sep_order = 'all' #'all', 'maj','accuy', 'cost' in separated models
    asc_order = 'desc'
    sep =  False # If a separated model
    sparse = False # If a sparse model

    for opt, arg in opts:
        if opt in ('-a'):
            asc_order = str(arg)
            assert asc_order in ('asc', 'desc'), 'option must be in (asc, desc)'
        elif opt in ('-h', '--help'):
            usage()
            sys.exit(0)
        elif opt in ('-l', '--lam', '--lambda'):
            lam = float(arg)
        elif opt in ('--maxsat', '--mxsat'):
            maxsat = True
        elif opt in ('-m', '--model'):
            ml = str(arg)
            assert ml in ('list', 'hybrid'), 'option must be in (list, hybrid)'
        elif opt in ('-n', '--node'):
            N = int(arg)
        elif opt in ('-o', '--order'):
            sep_order = str(arg)
            assert sep_order in ('all', 'maj', 'accuy', 'cost'), 'option must be in (all, maj, accuy, cost)'
        elif opt in ('--sep', '--separate'):
            sep = True
        elif opt in ('--spa', '--sparse'):
            sparse = True

        else:
            assert False, 'Unhandled option: {0} {1}'.format(opt, arg)


    return train_filename, N, lam, sep_order, asc_order, maxsat, sep, ml, sparse
#
#==============================================================================
def usage():
    """
        Prints help message.
    """

    print('Usage:', os.path.basename(sys.argv[0]), '[options] file')
    print('Options:')
    print('        -a                               The order)')
    print('                                         \'asc\': smallest/lowest first')
    print('                                         \'desc\': larggest/highest first')
    print('        -h, --help                       Print this usage message')
    print('        -l, --lam, --lambda==<float>     Value of lambda')
    print('        --maxsat, --mxsat                It is a MaxSAT model')
    print('        -m, --model==<string>            The model is a list, hybrid model')
    print('        -n, --nodes==<int>               The initial upper bound')
    print('        -o, --order==<string>            The class order in a separated model')
    print('                                         \'all\': output all possible orders')
    print('                                         \'maj\': output the order based on the majority of class')
    print('                                         \'accuy\': output the order based on the accuracy in each class')
    print('        --sep, --separate                It is a separated model')
    print('        --spa, --sparse                  It is a sparse model')

def compute(train_filename, num_nodes, lam, mis_weight=1, s_order='all', asc_order = 'desc', maxsat=True, sep=True, ml='list', sparse=False):

    # extract
    feature_dict, feature_names, feature_vars, data_features, data_classes, num_class = data_processing(train_filename)

    K = len(data_features[0])  # The number of features

    # To record the number of training items in each class
    class_count = [0 for i in range(num_class)]
    for i in range(len(data_classes)):
        for j in range(len(data_classes[i])):
            if data_classes[i][j] == 1:
                class_count[j] += 1
                continue

    # lambda value. Only used for sparse models:
    LAM = math.ceil(lam * len(data_features))


    if not sparse:
        # If not sparse, deal with inconsistent dataset
        feature_class_dict = {} # dictionary for counting classes of a feature set
        distinct_class = [] # the distinct classes

        # the distinct classes using one-hot encoding
        for c in range(num_class):
            distinct_class.append([0 if i != c else 1 for i in range(num_class)])

        # Counting classes of a feature set
        for i in range(len(data_features)):
            try:
                feature_class_dict[repr(data_features[i])][repr(data_classes[i])] += 1
            except:
                feature_class_dict[repr(data_features[i])] = {}

                for c in range(num_class):
                    feature_class_dict[repr(data_features[i])][repr(distinct_class[c])] = 0
                feature_class_dict[repr(data_features[i])][repr(data_classes[i])] += 1

        # Record the most class in a feature set
        for i in range(len(data_features)):
            max_class, max_class_count = [], 0
            for c in range(num_class):
                if feature_class_dict[repr(data_features[i])][repr(distinct_class[c])] > max_class_count:
                    max_class, max_class_count = list((distinct_class[c])), feature_class_dict[repr(data_features[i])][repr(distinct_class[c])]
            feature_class_dict[repr(data_features[i])]['max'] = max_class

        # Replace the class of a feature set if the original class of the feature set is not the most class of the feature set
        for i in range(len(data_features)):
            if list(data_classes[i]) != feature_class_dict[repr(data_features[i])]['max']:
                data_classes[i] = feature_class_dict[repr(data_features[i])]['max']
            assert list(data_classes[i]) == feature_class_dict[repr(data_features[i])]['max'], "Overlap! Index: {0}, {1}".format(i, j)

    # Compute the separated Model
    if sep:
        start_node = num_nodes # Assign a model size as the start point of computation
        sep_order = s_order # Order type: based on the number of items / accuracy / cost of each class

        # Sort data based on their class
        all_data_feature = [[] for i in range(num_class)]
        for i in range(len(data_features)):
            for j in range(num_class):
                if data_classes[i][j] == 1:
                    all_data_feature[j].append(list(data_features[i]))

        data = []
        for i in range(num_class):
            selected_features = []
            non_selected_features = []
            for j in range(num_class):
                if i != j:
                    non_selected_features += all_data_feature[j]
                else:
                    selected_features.append(all_data_feature[j])

            data.append(selected_features[0])

        data_dict = {class_ids + 1: data[class_ids][:] for class_ids in range(len(data))}

        # The class order for separated models
        all_class_order = []
        if ml == 'list':
            if sep_order == 'all' and num_class <= 3: # for all possible class orders
                for c in itertools.permutations(range(1, num_class + 1), num_class):
                    all_class_order.append(list(c))
            elif sep_order not in ('accuy', 'cost'):
                if sep_order == 'maj':
                    # sep_order == 'maj'
                    # Class order in the number of items in each class
                    # E.g. class 1 has 4 items, class 2 has 6 items, the order would be [6, 4]
                    sort_c_count = [sorted(class_count).index(x) for x in class_count]
                    for ii in range(len(class_count)):
                        new = sort_c_count[:]
                        new.pop(ii)
                        while sort_c_count[ii] in new:
                            sort_c_count[ii] += 1

                    maj_order = [None for j in range(len(class_count))]

                    for ii in range(len(class_count)):
                        maj_order[len(class_count) - 1 - sort_c_count[ii]] = ii + 1


                    if asc_order == 'desc': # class of more items first
                        all_class_order = [maj_order]
                    else:
                        assert asc_order == 'asc' # class of fewer items first
                        maj_order.reverse()
                        all_class_order = [maj_order]
                else:
                    assert sep_order == 'all' and num_class > 3
                    # Get 3 random orders
                    b_order = list(range(1, num_class + 1))
                    while True:
                        new_order = b_order[:]
                        random.shuffle(new_order)
                        if new_order not in all_class_order:
                            all_class_order.append(new_order)
                        if len(all_class_order) == 3:
                            break
                    assert len(all_class_order) == 3
        else:
            all_class_order.append(list(range(1, num_class+1))) # all classes

        all_order_detail = [] # Record the results of decision lists

        # The class order based on training accuracy
        # For step 1, the model goes through all classes respectively, and get the number of misclassification
        # The first decision list would be generated for the class with least misclassification
        # Repeat the step
        if sep_order in ('accuy', 'cost'):
            tem_order = list(range(1, num_class + 1)) # All possible classes for the first step
            for c_ids in range(num_class):
                all_class_order.append(tem_order[c_ids:] + tem_order[:c_ids])

            train_sol_rec = [{} for c_ids in range(num_class) ]
            train_sol_rec[0] = {c_ids + 1: {'mis_count': None, 'sol': None, 'best': False,
                                            'node': None, 'valid_node': None,'lit': None, 'time': 0,
                                            'data': None, 'cost': None}
                                for c_ids in range(num_class)}
            final_order = []

        done = False

        # For accuracy order only
        curr_pos = 0

        while not done:
            for class_order in all_class_order:
                if curr_pos == 0:
                    data = [data_dict[class_order[ii]][:] for ii in range(len(class_order))] # Data for the particular order
                elif sep_order in ('accuy', 'cost'):
                    data = [accu_new_data[class_order[ii]][:] for ii in range(len(class_order))]
                sols = [] # Record the solution for a particular order
                num_nodes = [] # Record the upper bound size
                num_valid_nodes = [] # Record the number of used nodes
                num_literal = [] # Number of literals

                # Record runtime
                all_runtime = [[] for i in range(num_class)]
                accu_runtime = 0

                # Upper bound for each class
                N = [start_node for i in range(num_class)]

                ubounds = [math.ceil((len(data_features[0]) + 1) * len(data[class_ids])) for class_ids in range(num_class)]

                # Go through all classes
                for class_ids in range(curr_pos, num_class):
                    found = False
                    solved = False

                    # compute
                    while not found:
                        DDL_solver = DLSolver(K, N[class_ids], data, class_ids=class_ids, LAM=LAM, mis_weight=mis_weight, maxsat=maxsat, sep=sep, ml=ml, sparse=sparse)
                        solutions = DDL_solver.solve()
                        all_runtime[class_ids].append(DDL_solver.runtime)
                        accu_runtime += DDL_solver.runtime

                        if solutions is not None:

                            for sol in solutions:

                                dl = DL(K, N[class_ids])
                                dl.parse_solution(sols=[sol])
                                dl.generate_node(cur_class=class_order[class_ids])

                                lit = dl.get_num_literals()
                                nof_miss = 0 #If not sparse


                                if sparse:
                                    # For upper bound
                                    nof_miss, misclas_ids = dl.get_miss(DDL_solver.M + DDL_solver.non_M, sol)
                                    nof_used = dl.valid_N
                                    cost = nof_used + int(math.ceil(nof_miss / LAM))
                                    ubounds[class_ids] = min(cost, ubounds[class_ids])
                                    if ubounds[class_ids] in (N[class_ids], nof_used):
                                        # either the formulation has reached the bound
                                        # or the actual solution did
                                        found = True
                                        solved = True
                                        sols.append(sol)

                                        filtered_data = DDL_solver.filtered_data[:]
                                        classified_data = DDL_solver.classified_data[:]

                                        if ml == 'list' and sparse: #Filter out classified data for the next step
                                            for ids in misclas_ids:
                                                for accu_num_ids in range(1, len(DDL_solver.num_data)):
                                                    if accu_num_ids == 1 and ids <= DDL_solver.num_data[accu_num_ids - 1]:
                                                        filtered_data[accu_num_ids - 1][ids - 1] = None
                                                        break
                                                    if ids <= DDL_solver.num_data[accu_num_ids] and ids > DDL_solver.num_data[accu_num_ids - 1]:
                                                        filtered_data[accu_num_ids][
                                                            ids - 1 - DDL_solver.num_data[accu_num_ids - 1]] = None
                                                        break

                                            for d in range(len(DDL_solver.filtered_data)):
                                                filtered_data[d] = list(filter(None, DDL_solver.filtered_data[d]))

                                            data = classified_data[:] + filtered_data[:]

                                            temp_accu_new_data = {}
                                            for ids in range(len(data)):
                                                temp_accu_new_data[class_order[ids]] = data[ids]


                                        num_nodes.append(N[class_ids])
                                        num_valid_nodes.append(dl.valid_N)
                                        num_literal.append(lit)

                                    else:
                                        if 10 < ubounds[class_ids] - N[class_ids]:
                                            if 10 < 2 * N[class_ids]:
                                                N[class_ids] += 10
                                            else:
                                                N[class_ids] *= 2
                                        else:
                                            N[class_ids] = ubounds[class_ids]


                                else:
                                    # If it not a sparse model, i.e. a perfect model
                                    found = True
                                    solved = True
                                    sols.append(sol)
                                    num_nodes.append(N[class_ids])
                                    num_valid_nodes.append(dl.valid_N)
                                    num_literal.append(lit)

                        else:
                            # If not solved, upper bound + 1
                            N[class_ids] += 1

                    # If the order is based on accuracy
                    # Record the result for each class respectively
                    # Only generate a decision list for the particular class in an order
                    if sep_order in ('accuy', 'cost'):
                        if not solved:
                            Detail = []
                            Detail.append(train_filename)
                            Detail.append(len(data_features))
                            Detail.append(len(data_features[0]))
                            Detail.append(num_class)
                            Detail.append('N')
                            Detail.append('not solved')
                            Detail.append('not solved')
                            Detail.append('not solved')
                            Detail.append('not solved')
                            Detail.append('not solved')
                            Detail.append('not solved')
                            Detail.append('not solved')
                            Detail.append('not solved')
                            Detail.append('not solved')

                        train_sol_rec[curr_pos][class_order[class_ids]] = {'mis_count': nof_miss, 'sol': sol, 'best': False,
                                                        'node': num_nodes[0], 'valid_node': dl.valid_N, 'lit': lit,
                                                        'time': sum(all_runtime[class_ids]), 'data': temp_accu_new_data,
                                                        'cost': nof_used + int(math.ceil(nof_miss / LAM)) if sparse else None}

                        break

                # If the order is not based on accuracy
                if sep_order not in ('accuy', 'cost'):
                    if solved:
                        dl = DL(K, sum(num_nodes))
                        dl.parse_solution(sols=sols,num_nodes=num_valid_nodes)
                        dl.generate_node(num_nodes=num_valid_nodes, class_order=class_order)

                        # Calculate the accuracy for training datasets
                        dl.mis_count = 0
                        dl.validate(data_features, data_classes, num_nodes=num_valid_nodes, num_class=num_class, ml=ml)
                        dl.generate_list(ml=ml, feature_dict=feature_dict, feature_names=feature_names)
                        train_accuracy = 1 - dl.mis_count / len(data_classes)

                        # Calculate the accuracy for testing datasets
                        dl.mis_count = 0

                        Detail = []
                        Detail.append(train_filename)
                        Detail.append(len(data_features))
                        Detail.append(len(data_features[0]))
                        Detail.append(num_class)
                        Detail.append('Y')
                        Detail.append(str(dl.valid_N) if dl.valid_N != 0 else 1)
                        Detail.append(sum(num_valid_nodes) - sum(num_literal) if sum(num_valid_nodes) != 0 else 1)
                        Detail.append(sum(num_literal))
                        Detail.append(round(train_accuracy * 100, 6))
                        #Detail.append(#round(test_accuracy * 100, 6))
                        Detail.append('%.6f' % accu_runtime)
                        Detail.append(class_order)
                        Detail.append(ml)
                        Detail.append(dl.node)

                        all_order_detail.append(Detail)
                    else:
                        Detail = []
                        Detail.append(train_filename)
                        Detail.append(len(data_features))
                        Detail.append(len(data_features[0]))
                        Detail.append(num_class)
                        Detail.append('N')
                        Detail.append('not solved')
                        Detail.append('not solved')
                        Detail.append('not solved')
                        Detail.append('not solved')
                        Detail.append('not solved')
                        Detail.append('not solved')
                        Detail.append(class_order)
                        Detail.append(ml)
                        Detail.append(dl.node)

                        all_order_detail.append(Detail)


            if sep_order not in ('accuy', 'cost'):
                done = True # Have gone through all orders
            else:
                # If the order is based on accuracy/cost
                if curr_pos == num_class-1: # Have gone through all classes
                    done = True
                    num_nodes = []
                    num_valid_nodes = []
                    lit = 0
                    accu_runtime = 0
                    class_order = final_order[:] + list(set(list(range(1, num_class + 1))) - set(final_order[:]))
                    final_order = class_order[:] # Final order based or accuracy
                    sols = []

                    # pick the matched class order
                    for c in range(num_class):
                        for cc in train_sol_rec[c].keys():
                            if c == num_class - 1:
                                train_sol_rec[c][cc]['best'] = True
                            accu_runtime += train_sol_rec[c][cc]['time']
                            if train_sol_rec[c][cc]['best']:
                                num_nodes.append(train_sol_rec[c][cc]['node'])
                                num_valid_nodes.append(train_sol_rec[c][cc]['valid_node'])
                                lit += train_sol_rec[c][cc]['lit']
                                sols.append(train_sol_rec[c][cc]['sol'])

                    #Calculate final accuracy
                    dl = DL(K, sum(num_nodes))
                    dl.parse_solution(sols=sols,num_nodes=num_valid_nodes)
                    dl.generate_node(num_nodes=num_valid_nodes, class_order=final_order)
                    dl.validate(data_features, data_classes, num_nodes=num_valid_nodes, num_class=num_class, ml=ml)

                    dl.generate_list(ml=ml, feature_dict=feature_dict, feature_names=feature_names)
                    train_accuracy = 1 - dl.mis_count / len(data_classes)

                    dl.mis_count = 0

                    # record the detail
                    Detail = []
                    Detail.append(train_filename)
                    Detail.append(len(data_features))
                    Detail.append(len(data_features[0]))
                    Detail.append(num_class)
                    Detail.append('Y')
                    Detail.append(sum(num_valid_nodes) if sum(num_valid_nodes) > 0 else 1)
                    Detail.append(sum(num_valid_nodes) - lit if sum(num_valid_nodes) - lit > 0 else 1) #'rule'
                    Detail.append(lit) # lit
                    Detail.append(round(train_accuracy * 100, 6))
                    Detail.append('%.6f' % accu_runtime)
                    Detail.append(final_order)
                    Detail.append(ml)
                    Detail.append(dl.node)

                    all_order_detail.append(Detail)

                else:
                    #If the final decision list has not been generated
                    # The order is based on accuracy
                    if sep_order == 'accuy':
                        # Record the class with least misclassification for a specific position is a decision list
                        if asc_order == 'desc':
                            min_mis_count = math.inf
                            min_class = None
                            for key in train_sol_rec[curr_pos].keys():
                                if train_sol_rec[curr_pos][key]['mis_count'] <= min_mis_count:
                                    min_mis_count = train_sol_rec[curr_pos][key]['mis_count']
                                    min_class = key

                            train_sol_rec[curr_pos][min_class]['best'] = True
                            best_class = min_class

                        else:
                            assert asc_order == 'asc'
                            # Record the class with most misclassification for a specific position is a decision list
                            max_mis_count = 0
                            max_class = None
                            for key in train_sol_rec[curr_pos].keys():
                                if train_sol_rec[curr_pos][key]['mis_count'] >= max_mis_count:
                                    max_mis_count = train_sol_rec[curr_pos][key]['mis_count']
                                    max_class = key

                            train_sol_rec[curr_pos][max_class]['best'] = True
                            best_class = max_class

                    else:
                        # The order is based on cost
                        assert sep_order == 'cost'
                        if asc_order == 'desc':
                            max_cost = 0
                            max_class = None
                            for key in train_sol_rec[curr_pos].keys():
                                if train_sol_rec[curr_pos][key]['cost'] >= max_cost:
                                    max_cost = train_sol_rec[curr_pos][key]['cost']
                                    max_class = key

                            train_sol_rec[curr_pos][max_class]['best'] = True
                            best_class = max_class
                        else:
                            assert asc_order == 'asc'
                            min_cost = math.inf
                            min_class = None
                            for key in train_sol_rec[curr_pos].keys():
                                if train_sol_rec[curr_pos][key]['cost'] <= min_cost:
                                    min_cost = train_sol_rec[curr_pos][key]['cost']
                                    min_class = key

                            train_sol_rec[curr_pos][min_class]['best'] = True
                            best_class = min_class


                    final_order.append(best_class)  # Order appends one best class
                    tem_order = list(range(1, num_class + 1)) # A possible order for the next step
                    n_order = list(set(tem_order) - set(final_order)) #The order for non-classified classes

                    # Possible orders for the other steps
                    all_class_order = []
                    for c_ids in range(len(n_order)):
                        all_class_order.append(final_order[:] + n_order[c_ids:] + n_order[:c_ids])


                    train_sol_rec[curr_pos+1] = {c_ids: {'mis_count': None, 'sol': None, 'best': False,
                                                    'node': None, 'valid_nodes': None, 'lit': None, 'time': 0,
                                                    'data': None, 'cost': None}
                                        for c_ids in n_order}

                    accu_new_data = train_sol_rec[curr_pos][best_class]['data'] #Used for new step
                    curr_pos += 1

        dl.node = all_order_detail[0][-1][:]
        dl.generate_list(feature_dict=feature_dict, feature_names=feature_names)
        print(dl.list)
    else:
        # For complete models
        #The upper bound
        N = [num_nodes]
        data = (data_features, data_classes) #training data
        sols = [] # Record the solutions
        num_nodes = [] # Number of unused and used nodes
        num_literal = [] # Number of literals excluding leaf nodes(rules)

        #Record the runtime
        all_runtime = []
        accu_runtime = 0

        found = False

        ubounds = [math.ceil((len(data_features[0]) + 1) * len(data_features))]

        # compute
        while not found:
            DDL_solver = DLSolver(K, N[0], data, LAM=LAM, mis_weight=mis_weight, maxsat=maxsat, sep=sep, ml=ml, sparse=sparse)

            solutions = DDL_solver.solve()
            all_runtime.append(DDL_solver.runtime)
            accu_runtime += DDL_solver.runtime

            if solutions is not None:

                for sol in solutions:

                    sols.append(sol)
                    num_nodes.append(N[0])

                    # Generate a decision list
                    dl = DL(K, N[0], C=num_class)
                    dl.parse_solution(sols=[sol])
                    dl.generate_node()

                    dl.validate(data_features, data_classes, num_nodes=[dl.valid_N])

                    nof_miss, misclas_ids = dl.get_miss(len(data[0]), sol)
                    nof_used = dl.valid_N
                    cost = nof_used + int(math.ceil(nof_miss / LAM))
                    ubounds[0] = min(cost, ubounds[0])
                    #assert nof_miss == dl.mis_count

                    # For size in sparse models
                    if ubounds[0] in (N[0], nof_used):
                        # either the formulation has reached the bound
                        # or the actual solution did
                        found = True
                        # Generate a decision lsit and calculate accuracy
                        dl.generate_list(feature_dict=feature_dict, feature_names=feature_names)
                        lit = dl.get_num_literals()
                        num_literal.append(lit)

                        dl.validate(data_features, data_classes, num_nodes=[dl.valid_N])
                        train_accuracy = 1 - dl.mis_count / len(data[0])

                        # For testing datasets
                        dl.mis_count = 0

                        # record detail
                        Detail = []
                        Detail.append(train_filename)
                        Detail.append(len(data_features))
                        Detail.append(len(data_features[0]))
                        Detail.append(num_class)
                        Detail.append('Y')
                        Detail.append(str(dl.valid_N))
                        Detail.append(dl.valid_N - lit if dl.valid_N != 1 else 1)
                        Detail.append(lit if dl.valid_N != 1 else 0)
                        Detail.append(round(train_accuracy * 100, 6))
                        Detail.append('%.6f' % accu_runtime)
                        Detail.append('Memory Peak')
                        Detail.append(ml)
                        Detail.append(dl.node)

                        print(dl.list)

                    else:
                        if 10 < ubounds[0] - N[0]:
                            if 10 < 2 * N[0]:
                                N[0] += 10
                            else:
                                N[0] *= 2
                        else:
                            N[0] = ubounds[0]

            else:
                N[0] += 1


if __name__ == "__main__":

        # get arguments from the terminal
        filename, N, lam, sep_order, asc_order, maxsat, sep, ml, sparse = parse_options()

        if sys.version_info.major == 2:
            sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

        # compute an optimal decision list with SAT
        compute(filename, N, lam, s_order=sep_order, asc_order = asc_order, maxsat=maxsat, sep=sep, ml=ml, sparse=sparse)

        if sys.version_info.major == 2:
            sys.stdout.close()
