#!/usr/bin/env python

import os

loc = ["check_0_universal", "check_1_beta", "check_2_fbeta"]
plots = ["jpsi_5_ATLAS", "jpsi_13_CMS", "psi2_13_CMS", "ups1_7_CMS", "ups2_13_CMS", "ups3_7_CMS"]

for i_l, l in enumerate(loc):
    for p in plots:
        os.system("cp %s/plots/%s_cs_pulls.pdf /home/mariana/Documents/2022_PhD_work/Phenom/1123_modelChecks/plots_model/%s_cs_pulls_v%d.pdf"%(l,p,p, i_l))
        os.system("cp %s/plots/%s_cs_devs.pdf /home/mariana/Documents/2022_PhD_work/Phenom/1123_modelChecks/plots_model/%s_cs_devs_v%d.pdf"%(l,p,p, i_l))

