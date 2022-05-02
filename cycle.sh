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
pdflatex plot.tex > res.txt
rm res.txt
cd ..

cp fit.txt fit_res/
cp -r plots fit_res/

# 2) call shell command to run fit over several values of a parameter
<<COMM
mkdir parscan
for file in l_*.txt
do
    echo "now doing file $file"
    cp $file param_list.txt	    
    root -l -q "c_test.C(1)" > log_root
    cd plots
    echo "doing plots"
    pdflatex plot.tex > res.txt
    rm res.txt
    cd ..
    echo "storing cosalpha"
    root -l -q "c_test_print.C"

    mkdir parscan/fit_res_$file
    mv fit.txt parscan/fit_res_$file/
    mv log_root parscan/fit_res_$file/
    mv cosa_scan.txt parscan/fit_res_$file/
    cp -r plots parscan/fit_res_$file/
done

echo "finished code"
COMM

# 3) call shell command to fit, get plots, then save everything in a folder
<<COMM
com="c_test_print.C"

root -l <<EOF
gSystem->Load("/home/mariana/local/lib/libLHAPDF.so");
gSystem->AddIncludePath("/home/mariana/local/include")
.x $com
.q
EOF

COMM
