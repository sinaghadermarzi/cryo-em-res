from __future__ import division

import csv
import sys
import copy
import operator
import numpy
import scipy

def get_amino_acid_data():
    csv_file_name = "amino_acid_data.csv"
    Data = {}
    with open(csv_file_name) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            name = row["ID"]
            del row["ID"]
            Data[name] = row
    return Data

def no_unique_chains(chains):
    seqs = []
    for ch in chains:
        seqs.append(ch["Sequence"])
    res = len(set(seqs))
    return res

def sum_chain_length(chains):
    sum = 0
    for row in chains:
        sum += int(row["Chain Length"])
    return sum

def average_chain_length(chains):
    sum = 0
    num = 0
    for row in chains:
        sum += int(row["Chain Length"])
        num +=1
    avg = sum / num
    return avg

def shortest_chain_length(chains):
    min = sys.maxint
    for row in chains:
        if (int(row["Chain Length"])<min):
            min = int(row["Chain Length"])
    return min

def longest_chain_length(chains):
    max = 0
    for row in chains:
        if (int(row["Chain Length"])>max):
            max = int(row["Chain Length"])
    return max

def num_chains(chains):
    return len(chains)

def stdev_chains(chains):
    ls = []
    for ch in chains:
        ls.append(len(ch["Sequence"]))
    sd = numpy.std(numpy.array(ls))
    return sd



def amino_acid_fractions_max(chains):
    num = {'A' : 0 ,'R' : 0 ,'N' : 0,'D' : 0 ,'C' : 0,'Q' : 0,'E' : 0,'G' : 0,'H': 0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0}
    row_num = []
    i = 0
    maxes = copy.deepcopy(num)
    for row in chains:
        sum = 0
        row_num = copy.deepcopy(num)
        seq = row["Sequence"]
        for ch in seq:
            row_num[ch] += 1
            sum += 1
        for k in num.keys():
            maxes[k] = max(maxes[k],row_num[k]/sum)
        i+=1
    return maxes

def amino_acid_fractions_min(chains):
    num = {'A' : 0 ,'R' : 0 ,'N' : 0,'D' : 0 ,'C' : 0,'Q' : 0,'E' : 0,'G' : 0,'H': 0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0}
    row_num = []
    i = 0
    mins = copy.deepcopy(num)
    for k in mins.keys():
        mins[k] = 1
    for row in chains:
        sum = 0
        row_num = copy.deepcopy(num)
        seq = row["Sequence"]
        for ch in seq:
            row_num[ch] += 1
            sum += 1
        for k in num.keys():
            mins[k] = min(mins[k],row_num[k]/sum)
        i+=1
    return mins

# def amino_acid_fractions_min(chains):
#     mm= sys.maxint
#     num = {'A' : mm,'R' : mm,'N' : mm,'D' : mm,'C' : mm,'Q' : mm,'E' : mm,'G' : mm,'H': mm,'I':mm,'L':mm,'K':mm,'M':mm,'F':mm,'P':mm,'S':mm,'T':mm,'W':mm,'Y':mm,'V':mm}
#     row_num = []
#     i = 0
#     for row in chains:
#         row_num[i] = copy(num)
#         seq = row["Sequence"]
#         for ch in seq:
#             row_num[i][ch] += 1
#         i+=1
#         for k in num.keys():
#             num[k] = min(num[k],row_num[i][k])
#     return num



def idx_felxiibility_sw_10_max(chains):
    am_data = get_amino_acid_data()
    min_cl = shortest_chain_length(chains)
    window_size = max(min_cl,10)
    as_max = 0
    for ch in chains:
        ch_max = 0
        temp = 0
        seq = ch['Sequence']
        for i in range(0, window_size):
            temp += float(am_data['Flexibility'][seq[i]])
        ch_max = temp
        for i in range (0,len(seq)-window_size):
            temp -= float(am_data['Flexibility'][seq[i]])
            temp += float(am_data['Flexibility'][seq[i+window_size]])
            if temp > ch_max:
                ch_max = temp
        if ch_max>as_max:
            as_max = ch_max
    return as_max

#########################################################################################



def max_avg_iupred_long_asamaksed_0_25(chains):
    ch_values =[]
    th_asa = 0.25
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_long']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))

        filtered_profbval = []
        asaquick_str = ch['ASAquick_Prediction']
        asaquick_spl_str = str.split(profbval_str,',')
        asaquick = []
        i = 0
        for s in asaquick_spl_str:
            asaquick.append(float(s))
        for i in range(0,len(profbval)):
            if asaquick[i]>th_asa:
                filtered_profbval.append(profbval[i])
        ch_values.append(numpy.mean(filtered_profbval))
    as_max = max(ch_values)
    return as_max


def max_iupred_long_sw_max_10_asamaksed_0_25(chains):
    th_asa = 0.25
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,10)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_long']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        filtered_profbval = []
        asaquick_str = ch['ASAquick_Prediction']
        asaquick_spl_str = str.split(profbval_str, ',')
        asaquick = []
        i = 0
        for s in asaquick_spl_str:
            asaquick.append(float(s))
        for i in range(0, window_size):
            if asaquick[i]>th_asa:
                temp += profbval[i]
        ch_max = temp
        for i in range(0,len(profbval)-window_size):
            if asaquick[i] > th_asa:
                temp -= profbval[i]
            if asaquick[i+window_size] > th_asa:
                temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
        ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max



def max_iupred_long_sw_max_20_asamaksed_0_25(chains):
    th_asa = 0.25
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_long']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        filtered_profbval = []
        asaquick_str = ch['ASAquick_Prediction']
        asaquick_spl_str = str.split(profbval_str, ',')
        asaquick = []
        i = 0
        for s in asaquick_spl_str:
            asaquick.append(float(s))
        for i in range(0, window_size):
            if asaquick[i]>th_asa:
                temp += profbval[i]
        ch_max = temp
        for i in range(0,len(profbval)-window_size):
            if asaquick[i] > th_asa:
                temp -= profbval[i]
            if asaquick[i+window_size] > th_asa:
                temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
        ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max


def max_iupred_long_sw_max_25_asamaksed_0_25(chains):
    th_asa = 0.25
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,25)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_long']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        filtered_profbval = []
        asaquick_str = ch['ASAquick_Prediction']
        asaquick_spl_str = str.split(profbval_str, ',')
        asaquick = []
        i = 0
        for s in asaquick_spl_str:
            asaquick.append(float(s))
        for i in range(0, window_size):
            if asaquick[i]>th_asa:
                temp += profbval[i]
        ch_max = temp
        for i in range(0,len(profbval)-window_size):
            if asaquick[i] > th_asa:
                temp -= profbval[i]
            if asaquick[i+window_size] > th_asa:
                temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
        ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max

#####################################################################################

def max_avg_iupred_short_asamaksed_0_25(chains):
    ch_values =[]
    th_asa = 0.25
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_short']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))

        filtered_profbval = []
        asaquick_str = ch['ASAquick_Prediction']
        asaquick_spl_str = str.split(profbval_str,',')
        asaquick = []
        i = 0
        for s in asaquick_spl_str:
            asaquick.append(float(s))
        for i in range(0,len(profbval)):
            if asaquick[i]>th_asa:
                filtered_profbval.append(profbval[i])
        ch_values.append(numpy.mean(filtered_profbval))
    as_max = max(ch_values)
    return as_max


def max_iupred_short_sw_max_10_asamaksed_0_25(chains):
    th_asa = 0.25
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,10)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_short']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        filtered_profbval = []
        asaquick_str = ch['ASAquick_Prediction']
        asaquick_spl_str = str.split(profbval_str, ',')
        asaquick = []
        i = 0
        for s in asaquick_spl_str:
            asaquick.append(float(s))
        for i in range(0, window_size):
            if asaquick[i]>th_asa:
                temp += profbval[i]
        ch_max = temp
        for i in range(0,len(profbval)-window_size):
            if asaquick[i] > th_asa:
                temp -= profbval[i]
            if asaquick[i+window_size] > th_asa:
                temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
        ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max



def max_iupred_short_sw_max_20_asamaksed_0_25(chains):
    th_asa = 0.25
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_short']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        filtered_profbval = []
        asaquick_str = ch['ASAquick_Prediction']
        asaquick_spl_str = str.split(profbval_str, ',')
        asaquick = []
        i = 0
        for s in asaquick_spl_str:
            asaquick.append(float(s))
        for i in range(0, window_size):
            if asaquick[i]>th_asa:
                temp += profbval[i]
        ch_max = temp
        for i in range(0,len(profbval)-window_size):
            if asaquick[i] > th_asa:
                temp -= profbval[i]
            if asaquick[i+window_size] > th_asa:
                temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
        ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max


def max_iupred_short_sw_max_25_asamaksed_0_25(chains):
    th_asa = 0.25
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,25)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_short']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        filtered_profbval = []
        asaquick_str = ch['ASAquick_Prediction']
        asaquick_spl_str = str.split(profbval_str, ',')
        asaquick = []
        i = 0
        for s in asaquick_spl_str:
            asaquick.append(float(s))
        for i in range(0, window_size):
            if asaquick[i]>th_asa:
                temp += profbval[i]
        ch_max = temp
        for i in range(0,len(profbval)-window_size):
            if asaquick[i] > th_asa:
                temp -= profbval[i]
            if asaquick[i+window_size] > th_asa:
                temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
        ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max





###########################################################################################





















###############################################################################

def max_avg_iupred_long(chains):
    ch_values =[]
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_long']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        ch_values.append(numpy.mean(profbval))
    as_max = max(ch_values)
    return as_max

def max_iupred_long_sw_max_10(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,10)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_long']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp
        for i in range (0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
            ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max

def max_iupred_long_sw_max_20(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_long']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp
        for i in range (0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
            ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max


def max_iupred_long_sw_max_25(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,25)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_long']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp
        for i in range (0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
            ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max


def avg_iupred_long_sw_max_10(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,10)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_long']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp
        for i in range (0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
            ch_values.append(ch_max)
    as_avg =numpy.mean(ch_values)
    return as_avg

def avg_iupred_long_sw_max_20(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_long']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp
        for i in range (0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
            ch_values.append(ch_max)
    as_avg =numpy.mean(ch_values)
    return as_avg


def avg_iupred_long_sw_max_25(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,25)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_long']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp
        for i in range (0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
            ch_values.append(ch_max)
    as_avg = numpy.mean(ch_values)
    return as_avg




def max_iupred_long_bin_frac_1(chains):
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        aa_values_str = ch['iupred_long']
        aa_values_spl_str = str.split(aa_values_str,',')
        aa_values = []
        i = 0
        for s in aa_values_spl_str:
            f = float(s)
            if f>0.5:
                aa_values.append(1)
            else:
                aa_values.append(0)
        frac= sum(aa_values)/len(aa_values)
        ch_values.append(frac)
    as_max = max(ch_values)
    return as_max

def max_iupred_long_bin_longest_seq1s(chains):
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        aa_values_str = ch['iupred_long']
        aa_values_spl_str = str.split(aa_values_str,',')
        aa_values = []
        i = 0
        for s in aa_values_spl_str:
            f = float(s)
            if f>0.5:
                aa_values.append(1)
            else:
                aa_values.append(0)
        count = 0
        max_consecutives = 0
        for a in aa_values:
            if a ==1:
                count+=1
                max_consecutives = max(count, max_consecutives)
            else:
                count = 0
        ch_values.append(max_consecutives)
    as_max = max(ch_values)
    return as_max




def avg_iupred_long_bin_frac_1(chains):
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        aa_values_str = ch['iupred_long']
        aa_values_spl_str = str.split(aa_values_str,',')
        aa_values = []
        i = 0
        for s in aa_values_spl_str:
            f = float(s)
            if f>0.5:
                aa_values.append(1)
            else:
                aa_values.append(0)
        frac= sum(aa_values)/len(aa_values)
        ch_values.append(frac)
    as_max = numpy.mean(ch_values)
    return as_max

def avg_iupred_long_bin_longest_seq1s(chains):
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        aa_values_str = ch['iupred_long']
        aa_values_spl_str = str.split(aa_values_str,',')
        aa_values = []
        i = 0
        for s in aa_values_spl_str:
            f = float(s)
            if f>0.5:
                aa_values.append(1)
            else:
                aa_values.append(0)
        count = 0
        max_consecutives = 0
        for a in aa_values:
            if a ==1:
                count+=1
                max_consecutives = max(count, max_consecutives)
            else:
                count = 0
        ch_values.append(max_consecutives)
    as_max = numpy.mean(ch_values)
    return as_max



#################################################################################


def max_avg_iupred_short(chains):
    ch_values =[]
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_short']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        ch_values.append(numpy.mean(profbval))
    as_max = max(ch_values)
    return as_max

def max_iupred_short_sw_max_10(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,10)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_short']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp
        for i in range (0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
            ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max

def max_iupred_short_sw_max_20(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_short']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp
        for i in range (0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
            ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max


def max_iupred_short_sw_max_25(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,25)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_short']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp
        for i in range (0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
            ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max


def avg_iupred_short_sw_max_10(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,10)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_short']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp
        for i in range (0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
            ch_values.append(ch_max)
    as_avg =numpy.mean(ch_values)
    return as_avg

def avg_iupred_short_sw_max_20(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_short']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp
        for i in range (0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
            ch_values.append(ch_max)
    as_avg =numpy.mean(ch_values)
    return as_avg


def avg_iupred_short_sw_max_25(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,25)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['iupred_short']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp
        for i in range (0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
            ch_values.append(ch_max)
    as_avg = numpy.mean(ch_values)
    return as_avg




def max_iupred_short_bin_frac_1(chains):
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        aa_values_str = ch['iupred_short']
        aa_values_spl_str = str.split(aa_values_str,',')
        aa_values = []
        i = 0
        for s in aa_values_spl_str:
            f = float(s)
            if f>0.5:
                aa_values.append(1)
            else:
                aa_values.append(0)
        frac= sum(aa_values)/len(aa_values)
        ch_values.append(frac)
    as_max = max(ch_values)
    return as_max

def max_iupred_short_bin_longest_seq1s(chains):
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        aa_values_str = ch['iupred_short']
        aa_values_spl_str = str.split(aa_values_str,',')
        aa_values = []
        i = 0
        for s in aa_values_spl_str:
            f = float(s)
            if f>0.5:
                aa_values.append(1)
            else:
                aa_values.append(0)
        count = 0
        max_consecutives = 0
        for a in aa_values:
            if a ==1:
                count+=1
                max_consecutives = max(count, max_consecutives)
            else:
                count = 0
        ch_values.append(max_consecutives)
    as_max = max(ch_values)
    return as_max




def avg_iupred_short_bin_frac_1(chains):
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        aa_values_str = ch['iupred_short']
        aa_values_spl_str = str.split(aa_values_str,',')
        aa_values = []
        i = 0
        for s in aa_values_spl_str:
            f = float(s)
            if f>0.5:
                aa_values.append(1)
            else:
                aa_values.append(0)
        frac= sum(aa_values)/len(aa_values)
        ch_values.append(frac)
    as_max = numpy.mean(ch_values)
    return as_max

def avg_iupred_short_bin_longest_seq1s(chains):
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        aa_values_str = ch['iupred_short']
        aa_values_spl_str = str.split(aa_values_str,',')
        aa_values = []
        i = 0
        for s in aa_values_spl_str:
            f = float(s)
            if f>0.5:
                aa_values.append(1)
            else:
                aa_values.append(0)
        count = 0
        max_consecutives = 0
        for a in aa_values:
            if a ==1:
                count+=1
                max_consecutives = max(count, max_consecutives)
            else:
                count = 0
        ch_values.append(max_consecutives)
    as_max = numpy.mean(ch_values)
    return as_max







#######################################

def max_avg_asaquick(chains):
    ch_values =[]
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['ASAquick_Prediction']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        ch_values.append(numpy.mean(profbval))
    as_max = max(ch_values)
    return as_max


































































































#################################################################################

def max_avg_profbval_asamaksed_0_25(chains):
    ch_values =[]
    th_asa = 0.25
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['ProfBVAL_Prediction']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))

        filtered_profbval = []
        asaquick_str = ch['ASAquick_Prediction']
        asaquick_spl_str = str.split(profbval_str,',')
        asaquick = []
        i = 0
        for s in asaquick_spl_str:
            asaquick.append(float(s))
        for i in range(0,len(profbval)):
            if asaquick[i]>th_asa:
                filtered_profbval.append(profbval[i])
        ch_values.append(numpy.mean(filtered_profbval))
    as_max = max(ch_values)
    return as_max


def max_profbval_sw_max_10_asamaksed_0_25(chains):
    th_asa = 0.25
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,10)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['ProfBVAL_Prediction']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        filtered_profbval = []
        asaquick_str = ch['ASAquick_Prediction']
        asaquick_spl_str = str.split(profbval_str, ',')
        asaquick = []
        i = 0
        for s in asaquick_spl_str:
            asaquick.append(float(s))
        for i in range(0, window_size):
            if asaquick[i]>th_asa:
                temp += profbval[i]
        ch_max = temp
        for i in range(0,len(profbval)-window_size):
            if asaquick[i] > th_asa:
                temp -= profbval[i]
            if asaquick[i+window_size] > th_asa:
                temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
        ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max



def max_profbval_sw_max_20_asamaksed_0_25(chains):
    th_asa = 0.25
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['ProfBVAL_Prediction']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        filtered_profbval = []
        asaquick_str = ch['ASAquick_Prediction']
        asaquick_spl_str = str.split(profbval_str, ',')
        asaquick = []
        i = 0
        for s in asaquick_spl_str:
            asaquick.append(float(s))
        for i in range(0, window_size):
            if asaquick[i]>th_asa:
                temp += profbval[i]
        ch_max = temp
        for i in range(0,len(profbval)-window_size):
            if asaquick[i] > th_asa:
                temp -= profbval[i]
            if asaquick[i+window_size] > th_asa:
                temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
        ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max


def max_profbval_sw_max_25_asamaksed_0_25(chains):
    th_asa = 0.25
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,25)
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['ProfBVAL_Prediction']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        filtered_profbval = []
        asaquick_str = ch['ASAquick_Prediction']
        asaquick_spl_str = str.split(profbval_str, ',')
        asaquick = []
        i = 0
        for s in asaquick_spl_str:
            asaquick.append(float(s))
        for i in range(0, window_size):
            if asaquick[i]>th_asa:
                temp += profbval[i]
        ch_max = temp
        for i in range(0,len(profbval)-window_size):
            if asaquick[i] > th_asa:
                temp -= profbval[i]
            if asaquick[i+window_size] > th_asa:
                temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
        ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max










###########################################################
def max_avg_prof_bval(chains):
    ch_values =[]
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['ProfBVAL_Prediction']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))



        ch_values.append(numpy.mean(profbval))
    as_max = max(ch_values)
    return as_max

def max_profbval_sw_max_10(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,10)
    as_max = 0
    ch_values = []
    for ch in chains:

        ch_max = 0
        temp = 0
        profbval_str = ch['ProfBVAL_Prediction']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max = temp

        for i in range(0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
        ch_values.append(ch_max)
        as_max = max(ch_values)
    return as_max

def max_profbval_sw_max_20(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    ch_values = []
    as_max = 0
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['ProfBVAL_Prediction']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max=temp
        for i in range(0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
        ch_values.append(ch_max)
    as_max = max(ch_values)
    return as_max


def max_profbval_sw_max_25(chains):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,25)
    as_max = 0
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['ProfBVAL_Prediction']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        for i in range(0, window_size):
            temp += profbval[i]
        ch_max=temp
        for i in range(0,len(profbval)-window_size):
            temp -= profbval[i]
            temp += profbval[i+window_size]
            if temp > ch_max:
                ch_max = temp
        if ch_max>as_max:
            as_max = ch_max
    return as_max





def max_avg_in_top10_percent_profbval(chains):
    perc = 30
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['ProfBVAL_Prediction']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        sorted_profbval = sorted(profbval)
        boundary = numpy.floor((1-perc/100)*len(profbval))
        sum =0
        count =0
        for i in range(int(boundary),len(profbval)):
            sum += sorted_profbval[i]
            count += 1
        avg = sum/count
    return avg




def min_avg_in_bottom10_precent_profbval(chains):
    perc = 30
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval_str = ch['ProfBVAL_Prediction']
        profbval_spl_str = str.split(profbval_str,',')
        profbval = []
        i = 0
        for s in profbval_spl_str:
            profbval.append(float(s))
        sorted_profbval = sorted(profbval)
        boundary = numpy.floor((perc/100)*len(profbval))
        sum =0
        count =0
        for i in range(0,int(boundary)):
            sum+=sorted_profbval[i]
            count+=1
        avg = sum/count
    return avg










def fraction_of_top3_flexibility_max(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['Flexibility'].items(), key=operator.itemgetter(1))
    #sorted_by_idx = list of 2-tuples of amino acids and their flexibility sorted by flexibility
    top_3 = sorted_by_idx[17:20]
    frac = amino_acid_fractions_max(chains)
    sum = 0
    for t in top_3:
        sum += frac[t[0]]
        amino_acid_fractions_max(chains)
    return sum

def fraction_of_top3_bvalue_max(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['B-Value'].items(), key=operator.itemgetter(1))
    #sorted_by_idx = list of 2-tuples of amino acids and their flexibility sorted by flexibility
    top_3 = sorted_by_idx[17:20]
    frac = amino_acid_fractions_max(chains)
    sum = 0
    for t in top_3:
        sum += frac[t[0]]
        amino_acid_fractions_max(chains)
    return sum


def fraction_of_top3_topidp_max(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['Top-IDP'].items(), key=operator.itemgetter(1))
    #sorted_by_idx = list of 2-tuples of amino acids and their flexibility sorted by flexibility
    top_3 = sorted_by_idx[17:20]
    frac = amino_acid_fractions_max(chains)
    sum = 0
    for t in top_3:
        sum += frac[t[0]]
        amino_acid_fractions_max(chains)
    return sum

def fraction_of_top3_fold_max(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['FoldUnfold'].items(), key=operator.itemgetter(1))
    #sorted_by_idx = list of 2-tuples of amino acids and their flexibility sorted by flexibility
    top_3 = sorted_by_idx[17:20]
    frac = amino_acid_fractions_max(chains)
    sum = 0
    for t in top_3:
        sum += frac[t[0]]
        amino_acid_fractions_max(chains)
    return sum


def fraction_of_top3_disprot_max(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['DisProt'].items(), key=operator.itemgetter(1))
    #sorted_by_idx = list of 2-tuples of amino acids and their flexibility sorted by flexibility
    top_3 = sorted_by_idx[17:20]
    frac = amino_acid_fractions_max(chains)
    sum = 0
    for t in top_3:
        sum += frac[t[0]]
        amino_acid_fractions_max(chains)
    return sum


####################################################


def fraction_of_bottom3_flexibility_max(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['Flexibility'].items(), key=operator.itemgetter(1))
    #sorted_by_idx = list of 2-tuples of amino acids and their flexibility sorted by flexibility
    top_3 = sorted_by_idx[0:2]
    frac = amino_acid_fractions_max(chains)
    sum = 0
    for t in top_3:
        sum += frac[t[0]]
        amino_acid_fractions_max(chains)
    return sum

def fraction_of_bottom3_bvalue_max(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['B-Value'].items(), key=operator.itemgetter(1))
    #sorted_by_idx = list of 2-tuples of amino acids and their flexibility sorted by flexibility
    top_3 = sorted_by_idx[0:2]
    frac = amino_acid_fractions_max(chains)
    sum = 0
    for t in top_3:
        sum += frac[t[0]]
        amino_acid_fractions_max(chains)
    return sum


def fraction_of_bottom3_topidp_max(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['Top-IDP'].items(), key=operator.itemgetter(1))
    #sorted_by_idx = list of 2-tuples of amino acids and their flexibility sorted by flexibility
    top_3 = sorted_by_idx[0:2]
    frac = amino_acid_fractions_max(chains)
    sum = 0
    for t in top_3:
        sum += frac[t[0]]
        amino_acid_fractions_max(chains)
    return sum

def fraction_of_bottom3_fold_max(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['FoldUnfold'].items(), key=operator.itemgetter(1))
    #sorted_by_idx = list of 2-tuples of amino acids and their flexibility sorted by flexibility
    top_3 = sorted_by_idx[0:2]
    frac = amino_acid_fractions_max(chains)
    sum = 0
    for t in top_3:
        sum += frac[t[0]]
        amino_acid_fractions_max(chains)
    return sum


def fraction_of_bottom3_disprot_max(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['DisProt'].items(), key=operator.itemgetter(1))
    #sorted_by_idx = list of 2-tuples of amino acids and their flexibility sorted by flexibility
    top_3 = sorted_by_idx[0:2]
    frac = amino_acid_fractions_max(chains)
    sum = 0
    for t in top_3:
        sum += frac[t[0]]
        amino_acid_fractions_max(chains)
    return sum


def max_in_chain_max_num_top3_flexibility_sw_20(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['Flexibility'].items(), key = operator.itemgetter(1))
    top_3 = sorted_by_idx[17:20]
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    as_max = 0
    considered_AAs = []
    for s in top_3:
        considered_AAs.append(s[0])
    for ch in chains:
        ch_max = 0
        temp = 0
        i = 0
        seq = ch['Sequence']
        for i in range(0, window_size):
            if seq[i] in considered_AAs:
                temp += 1
        ch_max=temp
        for i in range(0,len(seq)-window_size):
            if seq[i] in considered_AAs:
                temp -= 1
            if seq[i+window_size] in considered_AAs:
                temp += 1
            if temp > ch_max:
                ch_max = temp
        if ch_max>as_max:
            as_max = ch_max
    return as_max

def max_in_chain_max_num_top3_bvalue_sw_20(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['B-Value'].items(), key = operator.itemgetter(1))
    top_3 = sorted_by_idx[17:20]
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    as_max = 0
    considered_AAs = []
    for s in top_3:
        considered_AAs.append(s[0])
    for ch in chains:
        ch_max = 0
        temp = 0
        i = 0
        seq = ch['Sequence']
        for i in range(0, window_size):
            if seq[i] in considered_AAs:
                temp += 1
        ch_max=temp
        for i in range(0,len(seq)-window_size):
            if seq[i] in considered_AAs:
                temp -= 1
            if seq[i+window_size] in considered_AAs:
                temp += 1
            if temp > ch_max:
                ch_max = temp
        if ch_max>as_max:
            as_max = ch_max
    return as_max


def max_in_chain_max_num_top3_topidp_sw_20(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['Top-IDP'].items(), key = operator.itemgetter(1))
    top_3 = sorted_by_idx[17:20]
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    as_max = 0
    considered_AAs = []
    for s in top_3:
        considered_AAs.append(s[0])
    for ch in chains:
        ch_max = 0
        temp = 0
        i = 0
        seq = ch['Sequence']
        for i in range(0, window_size):
            if seq[i] in considered_AAs:
                temp += 1
        ch_max=temp
        for i in range(0,len(seq)-window_size):
            if seq[i] in considered_AAs:
                temp -= 1
            if seq[i+window_size] in considered_AAs:
                temp += 1
            if temp > ch_max:
                ch_max = temp
        if ch_max>as_max:
            as_max = ch_max
    return as_max

def max_in_chain_max_num_top3_fold_sw_20(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['FoldUnfold'].items(), key = operator.itemgetter(1))
    top_3 = sorted_by_idx[17:20]
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    as_max = 0
    considered_AAs = []
    for s in top_3:
        considered_AAs.append(s[0])
    for ch in chains:
        ch_max = 0
        temp = 0
        i = 0
        seq = ch['Sequence']
        for i in range(0, window_size):
            if seq[i] in considered_AAs:
                temp += 1
        ch_max=temp
        for i in range(0,len(seq)-window_size):
            if seq[i] in considered_AAs:
                temp -= 1
            if seq[i+window_size] in considered_AAs:
                temp += 1
            if temp > ch_max:
                ch_max = temp
        if ch_max>as_max:
            as_max = ch_max
    return as_max

def max_in_chain_max_num_top3_disprot_sw_20(chains):
    data = get_amino_acid_data()
    sorted_by_idx = sorted(data['DisProt'].items(), key = operator.itemgetter(1))
    top_3 = sorted_by_idx[17:20]
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,20)
    as_max = 0
    considered_AAs = []
    for s in top_3:
        considered_AAs.append(s[0])
    for ch in chains:
        ch_max = 0
        temp = 0
        i = 0
        seq = ch['Sequence']
        for i in range(0, window_size):
            if seq[i] in considered_AAs:
                temp += 1
        ch_max=temp
        for i in range(window_size,len(seq)-window_size):
            if seq[i] in considered_AAs:
                temp -= 1
            if seq[i+window_size] in considered_AAs:
                temp += 1
            if temp > ch_max:
                ch_max = temp
        if ch_max>as_max:
            as_max = ch_max
    return as_max

# def idx_bvalue_sw_min(chains):
#     am_data = get_amino_acid_data()
#     min_cl = shortest_chain_length(chains)
#     window_size  = max(min_cl,10)
#     sum = 0
#     seq = chains[0]["Sequence"]
#     for i in range(1,window_size):
#         sum+= float(am_data['Flexibility'][seq[i]])
#     for i in range(window_size)
#     return sum

# def max_avg_profbval(chains):
#     max_avg =0
#     for ch in chains:
#         profbval_str = ch['ProfBVAL_Prediction']
#         profbval_spl_str = str.split(profbval_str,',')
#         profbval = []
#         i = 0
#         for s in profbval_spl_str:
#             profbval.append(float(s))
#         avg = numpy.mean(profbval)
#         if avg > max_avg:
#             max_avg = avg
#     return max_avg


compute_feature = {

##TESTED FEATURES
    101:sum_chain_length,
    102:average_chain_length,
    103:shortest_chain_length,
    104:longest_chain_length,
##NOT TESTED FEATURES
    105:num_chains,
    106:stdev_chains,
    107:no_unique_chains,
    301:max_avg_iupred_long,
    302:max_iupred_long_sw_max_10,
    303:max_iupred_long_sw_max_20,
    304:max_iupred_long_sw_max_25,
    310:max_iupred_long_bin_frac_1,
    311:max_iupred_long_bin_longest_seq1s,
    320:avg_iupred_long_sw_max_10,
    321:avg_iupred_long_sw_max_20,
    322:avg_iupred_long_sw_max_25,
    323:avg_iupred_long_bin_frac_1,
    324:avg_iupred_long_bin_longest_seq1s,
    330:max_avg_iupred_short,
    331:max_iupred_short_sw_max_10,
    332:max_iupred_short_sw_max_20,
    333:max_iupred_short_sw_max_25,
    334:max_iupred_short_bin_frac_1,
    340:max_iupred_short_bin_longest_seq1s,
    350:avg_iupred_short_sw_max_10,
    351:avg_iupred_short_sw_max_20,
    352:avg_iupred_short_sw_max_25,
    353:avg_iupred_short_bin_frac_1,
    354:avg_iupred_short_bin_longest_seq1s,
    # 106:idx_felxiibility_sw_max,
    201:max_avg_prof_bval,
    202:max_profbval_sw_max_10,
    203:max_profbval_sw_max_20,
    204:max_profbval_sw_max_25,
    210:max_avg_in_top10_percent_profbval,
    211:min_avg_in_bottom10_precent_profbval,
    401:fraction_of_top3_flexibility_max,
    402:fraction_of_top3_bvalue_max,
    403:fraction_of_top3_topidp_max,
    404:fraction_of_top3_fold_max,
    405:fraction_of_top3_disprot_max,
    420:fraction_of_bottom3_flexibility_max,
    421:fraction_of_bottom3_bvalue_max,
    422:fraction_of_bottom3_topidp_max,
    423:fraction_of_bottom3_fold_max,
    424:fraction_of_bottom3_disprot_max,
    451:max_in_chain_max_num_top3_flexibility_sw_20,
    452:max_in_chain_max_num_top3_bvalue_sw_20,
    453:max_in_chain_max_num_top3_topidp_sw_20,
    454:max_in_chain_max_num_top3_fold_sw_20,
    455:max_in_chain_max_num_top3_disprot_sw_20,
    801:max_avg_asaquick,
    802:max_avg_profbval_asamaksed_0_25,
    803:max_profbval_sw_max_10_asamaksed_0_25,
    804:max_profbval_sw_max_20_asamaksed_0_25,
    805:max_profbval_sw_max_25_asamaksed_0_25,
    820:max_avg_iupred_long_asamaksed_0_25,
    821:max_iupred_long_sw_max_10_asamaksed_0_25,
    822:max_iupred_long_sw_max_20_asamaksed_0_25,
    823:max_iupred_long_sw_max_25_asamaksed_0_25,

}

