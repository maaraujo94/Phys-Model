#!/usr/bin/env python

import os

os.system("ls data_second_test > filelist.txt")

fin = open("filelist.txt")
fnames = fin.readlines()
fin.close()

for line in fnames:
    y = line.split('_')
    if y[0]=='LHCb' or y[2] == '7':
        os.system('cp data_second_test/%s data_third_test'%line[:-1])
    else:
        fin = open('data_second_test/%s'%line[:-1])
        fdata = fin.readlines()
        fin.close()
        fout = open('data_third_test/%s'%line[:-1],'w')
        for fl in fdata:
            y = fl.split()
            if float(y[3]) < 30:
                continue
            else:
                fout.write(fl)
        fout.close()
