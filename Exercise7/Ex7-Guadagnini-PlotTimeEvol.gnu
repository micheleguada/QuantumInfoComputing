# the maximum time is passed as command argument
# call this script as:
# gnuplot -c ../Ex7-Guadagnini-PlotTimeEvol.gnu tmax
#
tmax = ARG1

# set graphical options
set grid
set xrange [-5:5]
set palette defined (0 'blue',  1 'red')
set cbrange [0:tmax]
set cblabel "t [a.u.]"
set xlabel "x [a.u.]"; set ylabel "|{/Symbol Y}(x,t)|"
set arrow from 1,0 to 1,0.8 nohead lc rgb 'black' dt 2

# save the output
set term pdf color enhanced size 5,4
set output "Plot_TimeEvol.pdf"
plot "Tstep0000.dat" u 1:2 w l lt palette frac 1./12 title "",\
     "Tstep0100.dat" u 1:2 w l lt palette frac 2./12 title "",\
     "Tstep0200.dat" u 1:2 w l lt palette frac 3./12 title "",\
     "Tstep0300.dat" u 1:2 w l lt palette frac 4./12 title "",\
     "Tstep0400.dat" u 1:2 w l lt palette frac 5./12 title "",\
     "Tstep0500.dat" u 1:2 w l lt palette frac 6./12 title "",\
     "Tstep0600.dat" u 1:2 w l lt palette frac 7./12 title "",\
     "Tstep0700.dat" u 1:2 w l lt palette frac 8./12 title "",\
     "Tstep0800.dat" u 1:2 w l lt palette frac 9./12 title "",\
     "Tstep0900.dat" u 1:2 w l lt palette frac 10./12 title "",\
     "Tstep1000.dat" u 1:2 w l lt palette frac 11./12 title "",\
     "Tstep1024.dat" u 1:2 w l lt palette frac 12./12 title "",\
     

