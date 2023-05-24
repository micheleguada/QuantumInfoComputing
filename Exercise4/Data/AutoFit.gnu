datafile = ARG1
length = strlen(datafile)
method = substr(datafile, 0, length-4)
ff(x) = a*x+b
set key right bottom
set yrange [:2.2]
set xlabel "Log_{10}(N)"; set ylabel "Log_{10}(Time)"
#set title sprintf("%s%s", "Plot and Fit: ",datafile)
set grid
fit ff(x) datafile u (log10($1)):(log10($2)) via a,b
set label 1 "model: f(x) = a*x+b" at 1.7,1.5 font "arialbd,18"
set label 2 sprintf("  a = %3.4f", a) at 1.7,1.1 font "arialbd,18"
set label 3 sprintf("  b = %3.4f", b) at 1.7,0.7 font "arialbd,18"
set term pdf color enhanced size 4,4
set output sprintf("%s%s%s","Fit_",method,".pdf")
plot datafile u (log10($1)):(log10($2)) w p lc "red" pt 1 ps 1 title method, ff(x) w l lc "blue" lw 1.2 title sprintf("%s%s", method, " fit")
