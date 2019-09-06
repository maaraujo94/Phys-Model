#!/bin/sh

#call shell command to run fit and then plot and obtain a pdf with all plots

echo 'do fit before plotting? (1) for yes, (0) for no'
read doFit
com="c_test.C($doFit)"

root -l $com <<EOF
EOF

cd plots
echo 'now in plots directory'
pdflatex plot.tex > res.txt
rm res.txt
cd ..

evince plots/plot.pdf &
