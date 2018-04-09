# reads the result of evaluation of features and based on thresholds select a set of features
import csv
import numpy
import pandas
import scipy

from pandas import *
import numpy as np
#from libraries.settings import *
from scipy.stats.stats import pearsonr
import itertools

th_corr_between_features = 0.2
th_diff_rand = 0.2
th_p_value = 0.001
th_num_features = 20

features_file_name ='Train_Set_equal_dist.csv_features.csv_eval_res.csv'
feature_mean_diff = []
feature_p_value = []
feature_names = []

with open(features_file_name) as csvfile:
    reader = csv.DictReader(csvfile)
    field_names = reader.fieldnames
    curr_row = -1
    for row in reader:

        feature_mean_diff.append(abs(float(row['MeanX - MeanRand'])))
        feature_p_value.append(float(row['P-Value']))
        feature_names.append(row['Feature Name'])


#sorted_idx = sorted(range(len(feature_mean_diff)), key=lambda k: feature_mean_diff[k])

sorted_idx_rev = numpy.argsort(feature_mean_diff)
sorted_idx = sorted_idx_rev[::-1]
a = 7

df = pandas.read_csv('Train_Set_equal_dist.csv_features.csv')

# print(df)

# correlations = {}
# columns = df.columns.tolist()

# This line can be used to find the correlation between two features by their name
# scipy.stats.spearmanr(df.loc[:, col_a], df.loc[:, col_b])


selection = []
for i in sorted_idx:
#calculate the correlation between feature i and all in the selection
    if len(selection)>=th_num_features:
        break
    if (feature_mean_diff[i]>th_diff_rand) and (feature_p_value[i]<th_p_value):
        correlation_exist = 0
        for j in selection:
            rho, pval = scipy.stats.spearmanr(df.loc[:, feature_names[i]], df.loc[:, feature_names[j]])
            if rho > th_corr_between_features:
                correlation_exist = 1
                break
        if correlation_exist ==0:
            selection.append(i)


a=4
# for col_a, col_b in itertools.combinations(columns, 2):
#     correlations[col_a + '__' + col_b] = scipy.stats.spearmanr(df.loc[:, col_a], df.loc[:, col_b])

for f in selection:
    print feature_names[f], '\t', feature_mean_diff[f]




# result = DataFrame.from_dict(correlations, orient='index')
# result.columns = ['PCC', 'p-value']
#
# print(result.sort_index())


# for i in sorted_idx:
#     print feature_names[i],'\t',feature_mean_diff[i],'\n'