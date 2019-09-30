#!/bin/sh

#call shell command to run fit and then plot and obtain a pdf with all plots
<<COMM
echo "do fit before plotting? (1) for yes, (0) for no"
read doFit
com="c_test.C($doFit)"

root -l -q $com 

cd plots
echo "now in plots directory"
pdflatex plot.tex > res.txt
rm res.txt
cd ..

evince plots/plot.pdf &
COMM
#call shell command to run fit over several values of a parameter

#for file in b_*
for file in c_*.txt
do
    echo "now doing file $file"
    cp $file param_list.txt	    
    root -l -q "c_test.C(1)" > log_root
    cd plots
    echo "doing plots"
    pdflatex plot.tex > res.txt
    rm res.txt
    cd ..
    cp log_root parscan/log_$file
    cp fit.txt parscan/fit_$file
    cp -r plots parscan/plots_$file
done

echo "finished code"

