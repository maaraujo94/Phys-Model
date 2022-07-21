#!/usr/bin/env python

import os

min_v = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

fin = open('../model_final/fitModel.C')
fm = fin.readlines()
fin.close()

for i_c, ptv in enumerate(min_v):
    print i_c
    os.system('mkdir check_%d'%i_c)
    os.system('cp -r ../model_final/* check_%d'%i_c)
  
    fout = open('check_%d/fitModel.C'%i_c, 'w')
    for line in fm:
        if 'double ptmmin' in line:
            fout.write('double ptmmin = %.1f;\n'%ptv)
        else:
            fout.write(line)
    fout.close()
