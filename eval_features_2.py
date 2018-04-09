import random
import numpy
import scipy
import csv
from matplotlib import pyplot
from sklearn import linear_model as lm


def read_features(features_file_name):
    feature_data = []
    with open(features_file_name) as csvfile:
        reader = csv.DictReader(csvfile)
        field_names = reader.fieldnames
        curr_row = -1
        for row in reader:
            row_fl = {};
            for i in row.keys():
                if i != 'PDB ID':
                    row_fl[i] = float(row[i])
                else:
                    row_fl['PDB ID'] = row[i]
            feature_data.append(row_fl)
            curr_row += 1
        #        data.pop(0)
        return feature_data


def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(xrange(n), r))
    return tuple(pool[i] for i in indices)

#
# def get_feat_array(idx_list, feat_list):  # feat_list and idx_list are zero-based
#     data = []
#     mask = numpy.zeros(100)
#     for i in idx_list:
#         mask[i] = 1
#     with open(features_file_name) as csvfile:
#         reader = csv.reader(csvfile)
#         curr_row = -1
#         for row in reader:
#             if mask[curr_row] == 1:
#                 row_fl = map(float, [row[s] for s in feat_list])
#                 data.append(row_fl)
#             curr_row += 1
#         #        data.pop(0)
#         return numpy.array(data)


# def get_feat_array(idx_list, feat_list):  # feat_list and idx_list are zero-based
#     data = []
#     mask = numpy.zeros(100)
#     for i in idx_list:
#         mask[i] = 1
#     with open(features_file_name) as csvfile:
#         reader = csv.reader(csvfile)
#         curr_row = -1
#         for row in reader:
#             if mask[curr_row] == 1:
#                 row_fl = map(float, [row[s] for s in feat_list])
#                 data.append(row_fl)
#             curr_row += 1
#         #        data.pop(0)
#         return numpy.array(data)


def get_feature_subarray(feature_data, obj_list, feat):  # feat_list and idx_list are zero-based
    data = []
    num_objs = len(feature_data)
    mask = numpy.zeros(num_objs)
    for i in obj_list:
        mask[i] = 1
    curr_row = 0
    for row in feature_data:
        if mask[curr_row] == 1:
            data.append(row[feat])
        curr_row += 1
    #        data.pop(0)
    return numpy.array(data)



# n_samp = 1000
# samp_size = 40
# training_set_size = 70
#
# train_idx_rng = range(0, training_set_size)
# feature_idx_list = [8]
# rand_score_list = []
# X_score_list = []
#
# for i in range(1, n_samp + 1):
#     # generate a list of random indexes in between observations in X
#     rand_features = numpy.random.rand(training_set_size, 1)
#     idx_win = random_combination(train_idx_rng, samp_size)
#     X = get_feat_array(idx_win, feature_idx_list)
#     Y = get_feat_array(idx_win, [1])
#     X_rand = rand_features[idx_win, :]
#     (n_X_rows, n_X_cols) = X.shape
#     rho, pval = scipy.stats.spearmanr(X_rand, Y)
#     rand_score_list.append(rho)
#     rho, pval = scipy.stats.spearmanr(X, Y)
#     X_score_list.append(rho)
#
# pyplot.hist(rand_score_list, alpha=0.5, label='rand_score_list',bins=100)
# pyplot.hist(X_score_list, alpha=0.5, label='X_score_list',bins=100)
# pyplot.legend(loc='upper right')
# pyplot.savefig("a.png")
#
# print 'T test started'
# (t , p) = scipy.stats.ttest_rel(rand_score_list , X_score_list)
# print 'T_test finished'
#
# print ('T', t, '\n')
# print ('P', p, '\n')



#for each feature that we have read from the feature file
#evaluate feature by running T_Test and also print the histogram with filename feature_<feature_name>

#open csvfile

features_file_name = 'Train_Set_equal_dist.csv_features.csv'
feature_data = read_features(features_file_name)
num_training_objs = len(feature_data)
features_to_eval =[]
for f in feature_data[0].keys():
    if (f != 'Native_Resolution' and f != 'PDB ID'):
        features_to_eval.append(f)

rand_features = numpy.random.rand(num_training_objs, 1)
a = get_feature_subarray(feature_data, [2,3,4], 'feature_105')
n_samp = 1000
samp_size = 40

num_plot = 0
rand_sets = [[0]]
rand_sets_string = ['0']
rand_features = numpy.random.rand(num_training_objs, 1)
for i in range(1,n_samp+1):
    r = random_combination(range(0, num_training_objs), samp_size)
    rand_sets.append(r)
    rand_sets_string.append(','.join(str(e) for e in r))
# rand_sets = rand_sets+rand_sets
# rand_sets_string = rand_sets_string+rand_sets_string
unique_rand_sets = len(set(rand_sets_string))
if unique_rand_sets != len(rand_sets):
    print 'Error!\n'

with open(features_file_name+'_eval_res.csv','wb') as out_csv:
    field_names = ['Feature Name', 'MeanX - MeanRand' , 'P-Value']
    writer = csv.DictWriter(out_csv, field_names)
    writer.writeheader()
    for f in features_to_eval:
        out_row = {}
        rand_score_list = []
        X_score_list = []
        # evaluate feature
        for i in range(1, n_samp + 1):
            # generate a list of random indexes in between observations in X

            idx_win = rand_sets[i]
            X = get_feature_subarray(feature_data,idx_win, f)
            Y = get_feature_subarray(feature_data,idx_win, 'Native_Resolution')
            X_rand = rand_features[idx_win, :]
            # (n_X_rows, n_X_cols) = X.shape
            rho, pval = scipy.stats.spearmanr(X_rand, Y)
            rand_score_list.append(rho)
            rho, pval = scipy.stats.spearmanr(X, Y)
            X_score_list.append(rho)
        # create a row with feature name , t_stat, p_value
        samp_idx = random_combination(range(0, n_samp), 100)
        X_ar = numpy.array(X_score_list)
        rand_ar = numpy.array(rand_score_list)
        X_score_list_samp = X_ar[list(samp_idx)]
        rand_score_list_samp = rand_ar[list(samp_idx)]
        (t, p) = scipy.stats.ttest_rel(rand_score_list_samp, X_score_list_samp)
        diff=numpy.mean(rand_score_list_samp)-numpy.mean(X_score_list_samp)
        out_row['Feature Name'] = f
        out_row['MeanX - MeanRand'] = diff
        out_row['P-Value'] = p
        writer.writerow(out_row)
        # save histogram in a file with the feature_<feature_name> as name
        pyplot.figure(num_plot)
        pyplot.hist(rand_score_list_samp,range=[-0.5, 0.5], alpha=0.5, label='rand_score_list', bins=25)
        pyplot.hist(X_score_list_samp,range=[-0.5, 0.5], alpha=0.5, label='X_score_list', bins=25)
        pyplot.legend(loc='upper right')
        pyplot.savefig("feature_"+str(f)+".png")
        num_plot += 1
