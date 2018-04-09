import openpyxl as opx
workbook = opx.load_workbook('Prepared_Data.xlsx')
#sh_assemblies = workbook['Assemblies(2015-2017)']
sh_train_chains = workbook['Train_Chains']


for i in range (2,765):
    # Open the file with read only permit
    fileaddress = "ASA_Quick_Res\\asaq."+str(i)+".fasta\\asaq.pred"
    f = open(fileaddress)
    # use readline() to read the first line 

    # use the read line to read further.
    # If the file is not empty keep reading one line
    # at a time, till the file is empty

    res_list = list()
    line = f.readline()
    while line:
        line = line.rstrip('\n')
        line_cols = line.split(" ")
        res_list.append(line_cols[3])
        # use realine() to read next line
        line = f.readline()
    f.close()
    res_list = res_list[:-1]
    entry_res_string = ",".join(res_list)
    sh_train_chains.cell(row=i,column=10).value = entry_res_string

workbook.template = False
workbook.save('Prepared_Data.xlsx')