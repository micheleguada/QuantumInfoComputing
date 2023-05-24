# gnuplot script for timings

# graphical options
#set logscale xy
set xlabel "N (# of spins)"; set ylabel "log_2(time [s])"      
set xrange [1:13]
set key top left
set grid
set format y "%.3g"
set xtics (2,3,4,5,6,7,8,9,10,11,12)
set ytics (-10,-8,-6,-4,-2,0,2,4,6,8,10)
set yrange [-10:10]

# fit 
f(x) = a + b*x; a = 5; b = 0.1;
fit [8:12] f(x) "Timings.txt" u 1:(log($2)/log(2)) via a,b

# set box with fitted parameters
set style textbox opaque noborder
set obj 10 rect from graph 0.4,graph 0.95 to graph 0.6,graph 0.70
set label 1 "f(N) = a + b*N"             at graph 0.42,0.9 boxed font "arialbd,14"
set label 2 sprintf("a = %3.4g", a)      at graph 0.45,0.86 boxed font "arialbd,14"
set label 3 sprintf("b = %3.4g", b)      at graph 0.45,0.82 boxed font "arialbd,14"
set label 4 "fit range: [8:12]"          at graph 0.42,0.75 boxed font "arialbd,14"

# save the picture
set term pdf color enhanced size 6,4
set output "Timings.pdf"
plot "Timings.txt" u 1:(log($2)/log(2)) w p lc "red" pt 7 title "Timings", f(x) w l lc "blue" title "linear fit" 
