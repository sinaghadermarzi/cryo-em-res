

ofile = open("dump.csv","w")
for i in range (2,765):
    # Open the file with read only permit
    fileaddress = "peripheral_codes\\ASA_Quick_Res\\asaq."+str(i)+".fasta\\rasaq.pred"
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
        value = (float(line_cols[2])/2)+0.5
        if value<0:
            print 'ERROR!'
        res_list.append(str(value))
        # res_list.append(line_cols[2])
        # use realine() to read next line
        line = f.readline()
    res_list = res_list[:-1]
    entry_res_string = ",".join(res_list) + '\n'
    ofile.write(entry_res_string)
    f.close()
ofile.close

