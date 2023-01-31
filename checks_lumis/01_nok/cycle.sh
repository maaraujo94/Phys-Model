#!/bin/sh

# currently two different methods available to run fit
# 1) fit, get plots, save everything in a folder "fit_res"
# 2) using [par]_[val].txt files, fit over a cycle of values and save results
# 3) run just once to store x values (c_test_print.C)

# 1) call shell command to fit, get plots, then save everything in a folder

echo "doing fitting"
com="c_test.C(0)"

root -l <<EOF
gSystem->Load("/home/mariana/local/lib/libLHAPDF.so");
.include /home/mariana/local/include
.x $com
.q
EOF

cd plots
echo "now in plots directory"
pdflatex plot_dev.tex > res.txt
pdflatex plot_pull.tex > res.txt
rm res.txt
cd ..
