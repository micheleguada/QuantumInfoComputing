##### Plot the difference between the first 80 eigenvalues #####

## Set plot options
set key right top
set yrange [-0.6:0.2]
#set xrange [0:1]
set xlabel "Exact eigenvalue"; set ylabel "Relative Error (%)"
set grid
#set logscale x
#set logscale y

## Save the plot
set term pdf color enhanced size 6,4
set output "EigenvaluesError.pdf"
plot "Eigenvalues.dat" u 1:((($1-$2)/$1)*100) w lp lc "red" pt 7 ps 0.5 title "Relative Error (%)", 0 w l lc "blue" title ""
