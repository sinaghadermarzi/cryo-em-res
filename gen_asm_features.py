from feature_implement import *
import csv
import sys

#this program reads an input csv file containing PDB ID of proteins that we want to compute features for
#then using the existing feature file that contains the chain features it creates a per assembly feature file
chain_data_file_name = 'Chain_Data.csv'
working_data_file_name = 'Train_Set_equal_dist.csv'
compute_all_features = 1

# features_to_compute = [
#     101,
#     102,
#     103,
#     104,
#     105
# ]






def find_chains(PDBID):
    found_chains = []
    with open(chain_data_file_name) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["PDB ID"] == PDBID:
                found_chains.append(row)
    return found_chains





#for each element in working data look for the associated id in the chain data and bring all data into a dictionary and based on that compute the desired feature
#feature__column = 3
if compute_all_features:
    features_to_compute = compute_feature
with open(working_data_file_name) as csvfile, open(working_data_file_name+'_features.csv','wb') as out_csv:
    reader = csv.DictReader(csvfile)
    field_names = reader.fieldnames
    for f in features_to_compute:
        field_names = field_names + ["feature_" + str(f)]
    for f in range(901,941):
        field_names = field_names + ["feature_" + str(f)]
    writer = csv.DictWriter(out_csv, field_names)
    writer.writeheader()
    ams = get_amino_acid_data()['Flexibility'].keys()
    for row in reader:
        PDBID = row["PDB ID"]
        chains = find_chains(PDBID)
        out_row = row
        for f in features_to_compute:
            feature = compute_feature[f](chains)
            out_row.update({"feature_"+str(f): str(feature)})
        first_am_feature = 900
        num = amino_acid_fractions_max(chains)
        i = 901
        for am in num:
            feature = num[am]
            out_row["feature_"+str(i)]= str(feature)
            i+=1
        num = amino_acid_fractions_min(chains)
        for am in num:
            feature = num[am]
            out_row["feature_"+str(i)]= str(feature)
            i+=1
        writer.writerow(out_row)
        #create new dictionary with row and feature and save it to the output file




