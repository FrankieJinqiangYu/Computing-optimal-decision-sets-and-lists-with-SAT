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
import gzip


#
#==============================================================================
def data_processing(file_name):
    """

        Function to process data with binary features

    """
    data = gzip.open(file_name, 'rt')

    feature_names = data.readline().strip("\n").split(",")
    feature_vars = [[] for i in range(len(feature_names))]
    num_examples = 0

    for i in range(len(feature_names)):
        feature_names[i] = feature_names[i].strip("\"").replace('"', '')

    end = False
    while not end:
        line = data.readline().strip("\n").split(",")
        if line == ['']:
            end = True
            data.close()
        else:
            for i, f in enumerate(line):
                if f not in feature_vars[i]:
                    feature_vars[i].append(f)
            num_examples += 1

    # sort feature
    for i in range(len(feature_vars)):
        feature_vars[i].sort()

    num_features = 0
    num_class = len(feature_vars[-1])
    for i in range(len(feature_names) - 1):
        c = len(feature_vars[i]) if len(feature_vars[i]) > 2 else 1
        num_features += c

    data_features = [ [0 for j in range(num_features)] for i in range(num_examples)]
    data_classes =[ [0 for j in range(num_class)] for i in range(num_examples)]

    data = gzip.open(file_name, 'rt')
    data.readline()  # Skip the header
    end = False
    curr_exmaple_index = 0

    feature_dict = {} # for record the exact values in a feature (not 1 or 0)

    while not end:
        line = data.readline().strip("\n").split(",")
        if line == ['']:
            end = True
            data.close()
        else:
            for i, f in enumerate(line[:-1]):
                num_prev_vars = 0
                for j in range(i):
                    num_prev_vars += len(feature_vars[j]) if len(feature_vars[j]) > 2 else 1

                if len(feature_vars[i]) > 2: # the current feature has 3 or more distinct values
                    curr_f_index = feature_vars[i].index(f) + num_prev_vars
                    data_features[curr_exmaple_index, curr_f_index] = 1
                else:
                    curr_f_index = num_prev_vars
                    if curr_exmaple_index == 0:
                        first_line = line
                        data_features[0][curr_f_index] = 0 if str(first_line[i]).strip().upper() in ('FALSE', '0', 'NO', 'WEAK', 'NORMAL') else 1

                        try:
                            curr_f_exact_value = int(first_line[i].strip())
                        except:
                            curr_f_exact_value = str(first_line[i]).strip()

                        feature_dict[feature_names[i]] = {}
                        feature_dict[feature_names[i]][data_features[0][curr_f_index]] = curr_f_exact_value

                    else:
                        data_features[curr_exmaple_index][curr_f_index] = data_features[0][curr_f_index] if line[i] == first_line[i] else (data_features[0][curr_f_index] + 1) % 2

                        try:
                            curr_f_exact_value = int(line[i].strip())
                        except:
                            curr_f_exact_value = str(line[i]).strip()

                        if len(feature_dict[feature_names[i]]) == 1 and line[i] != first_line[i]:
                            feature_dict[feature_names[i]][data_features[curr_exmaple_index][curr_f_index]] = curr_f_exact_value

            if feature_names[-1] not in feature_dict:
                feature_dict[feature_names[-1]] = {}

            try:
                curr_f_exact_value = int(line[-1].strip())
            except:
                curr_f_exact_value = str(line[-1]).strip()

            feature_dict[feature_names[-1]][feature_vars[-1].index(line[-1])] = curr_f_exact_value

            data_classes[curr_exmaple_index][feature_vars[-1].index(line[-1])] = 1
            curr_exmaple_index += 1


    return feature_dict, feature_names, feature_vars, data_features, data_classes, num_class