optflag = ARG1
set xlabel "Matrix size N"
set ylabel "Exec. time [s]"
#set title sprintf("%s%s%s", "Different methods comparison (-",optflag," flag)")
set grid
set xrange [40:2500]
set yrange [:200]
set logscale yx
set format y "%.3g"
set key bottom right notitle nobox
set term pdf color size 4,4
set output sprintf("%s%s%s","Methods",optflag,".pdf")
plot sprintf("%s%s%s","1stLoop",optflag,".txt") w lp title "1stLoop", sprintf("%s%s%s","2ndLoop",optflag,".txt") w lp title "2ndLoop", sprintf("%s%s%s","MATMUL",optflag,".txt") w lp title "MATMUL"
