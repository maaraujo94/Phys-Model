#!/usr/bin/env python

import os

pdfset = ["CT18NNLO", "NNPDF31_nnlo_as_0118", "MMHT2014nnlo68cl", "ABMP16_5_nnlo"]

fin = open('../model_final/aux_func.C')
fm = fin.readlines()
fin.close()

for i_c, pdf in enumerate(pdfset):
    print i_c
    os.system('mkdir check_%d'%i_c)
    os.system('cp -r ../model_final/* check_%d'%i_c)
  
    fout = open('check_%d/aux_func.C'%i_c, 'w')
    for line in fm:
        if 'PDF *pdf_ct' in line:
            fout.write('PDF *pdf_ct = mkPDF("%s", 0);\n'%pdf)
        else:
            fout.write(line)
    fout.close()
