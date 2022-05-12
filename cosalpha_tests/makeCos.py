#!/usr/bin/env python

import os

rho_v = [2,3]
delta_v = [0,1]
n_type = 3

fin = open('../model_beta/param_list.txt')
params = fin.readlines()
fin.close()

i_c = 0
for t in range(n_type):
    for rho in rho_v:
        for delta in delta_v:
            os.system('mkdir check_%d'%i_c)
            os.system('cp ../model_beta/*.C check_%d'%i_c)
            os.system('cp ../model_beta/cycle.sh check_%d'%i_c)
            os.system('cp -r ../model_beta/plots check_%d'%i_c)
  
            fout = open('check_%d/param_list.txt'%i_c, 'w')
            for pl in params:
                if 'rho' in pl:
                    fout.write('shape rho %d 0.5 1\n'%rho)
                elif 'delta' in pl:
                    fout.write('shape delta %d 0.5 1\n'%delta)
                else:
                    fout.write(pl)
            fout.close()

            os.system('cp aux_func_%d.C check_%d/aux_func.C'%(t, i_c))
            
            i_c += 1
