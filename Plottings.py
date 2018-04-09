from __future__ import division

import csv
import sys
import copy
import operator
import numpy
import scipy
from feature_implement import *
from matplotlib import pyplot

# am_data = get_amino_acid_data()
# flex = am_data['FoldUnfold'].values()
# flex_f = []
# for i in flex:
#     flex_f.append(float(i))
# pyplot.figure()
# pyplot.hist(flex_f, alpha=0.5, label='flex', bins=10)
# pyplot.show()


chain_data_file_name = 'Chain_Data.csv'
working_data_file_name = 'Train_Set_equal_dist.csv'
def find_chains(PDBID):
    found_chains = []
    with open(chain_data_file_name) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["PDB ID"] == PDBID:
                found_chains.append(row)
    return found_chains


with open(working_data_file_name) as csvfile, open(working_data_file_name+'_features.csv','wb') as out_csv:
    reader = csv.DictReader(csvfile)
    good_res_profs= []
    med_res_profs = []
    bad_res_profs = []
    for row in reader:
        PDBID = row["PDB ID"]
        chains = find_chains(PDBID)
        out_row = row
        res =float(row['Native_Resolution'])
        if res<4.45:
            for ch in chains:
                profbval_str = ch['ProfBVAL_Prediction']
                profbval_spl_str = str.split(profbval_str, ',')
                profbval = []
                i = 0
                for s in profbval_spl_str:
                    profbval.append(float(s))
                good_res_profs.append(profbval)
        if res > 4.45 and res<9.85:
            for ch in chains:
                profbval_str = ch['ProfBVAL_Prediction']
                profbval_spl_str = str.split(profbval_str, ',')
                profbval = []
                i = 0
                for s in profbval_spl_str:
                    profbval.append(float(s))
                bad_res_profs.append(profbval)


# count = 0
# pyplot.figure(1)
# for i in good_res_profs:
#     if count%70 ==0:
#         t = numpy.array(range(0,len(i)))/len(i)
#         s = numpy.array(i)
#         # t = numpy.arange(0.0, 5.0, 0.01)
#         # s = np.cos(2 * np.pi * t)
#         line, = pyplot.plot(t, s, lw=1,color='C2')
#     count+=1
# pyplot.show(block=False)


pyplot.figure(1)
k = 0
count = 0
for i in bad_res_profs:
    if k%20 ==0:
        t = numpy.array(range(0,len(i)))
        s = numpy.array(i)
        # t = numpy.arange(0.0, 5.0, 0.01)
        # s = np.cos(2 * np.pi * t)
        line, = pyplot.plot(t, s, lw=1,color='tab:blue', alpha=0.5)

        t = numpy.array(range(0,len(good_res_profs[k])))
        s = numpy.array(good_res_profs[k])
        line, = pyplot.plot(t, s, lw=1 ,color='tab:red',alpha=0.5)
        count+=1
    k+=1
pyplot.show(block=False)
print 'count = ', count
# input('press')
pyplot.figure(3)
pyplot.show()

