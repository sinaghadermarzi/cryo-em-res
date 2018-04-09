import openpyxl as opx
import re
workbook = opx.load_workbook('Prepared_Data.xlsx')
#sh_assemblies = workbook['Assemblies(2015-2017)']
sh_train_chains = workbook['Train_Chains']


for i in range(2,356)+range(357,765):
    # Open the file with read only permit
    fileaddress = "ProfBVAL_Results\\"+str(i)+".bVal"
    f = open(fileaddress)
    # use readline() to read the first line 

    # use the read line to read further.
    # If the file is not empty keep reading one line
    # at a time, till the file is empty

    res_list = list()
    line = f.readline()
    line = f.readline()
    while line:
        line = line.rstrip('\n')
        line_cols = line.split("\t")      
        res_list.append(re.sub(' ','',line_cols[3]))
        # use realine() to read next line
        line = f.readline()
    f.close()
#    res_list = res_list[:-1]
    entry_res_string = ",".join(res_list)
    sh_train_chains.cell(row=i,column=11).value = entry_res_string

workbook.template = False
workbook.save('Prepared_Data.xlsx')
