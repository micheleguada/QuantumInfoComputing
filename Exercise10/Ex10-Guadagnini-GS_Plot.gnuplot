# gnuplot script to plot the ground state eigenvalue obtained with 
# RSRG method in function of Lambda for different N. 
# As a comparison, also the mean field solution is added.
#
# To run: gnuplot Ex10-Guadagnini-GS_Plot.gnuplot 

# graphical options
set xlabel "{/Symbol l}"; set ylabel "e_{GS}({/Symbol l}) "
set key bottom left
set grid
# set format y "%.3g"

# files selection
FILES   = system("ls -1 Data/EigVals_*.dat")
NLABELS = system("ls -1 Data/EigVals_*.dat | cut -c16-16")
TLABELS = system("ls -1 Data/EigVals_*.dat | cut -c20-20")

## plot all in one: ##

# line styles
set style line  1 dt 1 lw 1 lc "black"
set style line  2 dt 2 lw 1 lc rgb '#440154' # dark purple
set style line  3 dt 4 lw 1 lc rgb '#472c7a' # purple
set style line  4 dt 5 lw 1 lc rgb '#3b518b' # blue
set style line  5 dt 2 lw 1 lc rgb '#2c718e' # blue
set style line  6 dt 4 lw 1 lc rgb '#21908d' # blue-green
set style line  7 dt 5 lw 1 lc rgb '#27ad81' # green
set style line  8 dt 2 lw 1 lc rgb '#5cc863' # green
set style line  9 dt 4 lw 1 lc rgb '#aadc32' # lime green
set style line 10 dt 5 lw 1 lc rgb '#fde725' # yellow
set style line 11 dt 1 lw 1 lc "red"

outname = "Plots/GroundStatesAll.pdf"
set term pdf color enhanced size 5,4
set output outname
plot "Data/Additional_EigVals_N02_T02.dat" u 3:4 w l ls 1 title sprintf( "N=%s; T=%s", "2", "2" ) ,\
     for [idx = 1:9] word(FILES,idx) u 3:4 w l ls (idx+1) title sprintf( "N=%s; T=%s", word(NLABELS,idx), word(TLABELS,idx) ) ,\
     "Data/MeanFieldGS.dat" u 1:2 w l ls 11 title " MeanField"
     
outname = "Plots/GroundStatesAll_Zoomed.pdf"
set xrange [0.4:1.4]
set term pdf color enhanced size 5,4
set output outname
plot "Data/Additional_EigVals_N02_T02.dat" u 3:4 w l ls 1 title sprintf( "N=%s; T=%s", "2", "2" ) ,\
     for [idx = 1:9] word(FILES,idx) u 3:4 w l ls (idx+1) title sprintf( "N=%s; T=%s", word(NLABELS,idx), word(TLABELS,idx) ) ,\
     "Data/MeanFieldGS.dat" u 1:2 w l ls 11 title " MeanField"

unset xrange     
## other plots: ##

# line styles
set style line  1 dt 1 lw 2 lc "blue"
set style line  2 dt 2 lw 2 lc "web-green"
set style line  3 dt 4 lw 2 lc "red"

set style line 11 dt 1 lw 2 lc "black"

# plot with same N
do for [jj = 1:9:3] {

    outname = sprintf("Plots/GroundStates_N%s.pdf", word(NLABELS,jj))  # Note: check that the folder "Plots" exists before running
    
    # save the picture
    set term pdf color enhanced size 5,4
    set output outname
    plot for [idx = 1:3] word(FILES,idx-1+jj) u 3:4 w l ls idx title sprintf( "N=%s; T=%s", word(NLABELS,jj), word(TLABELS,idx) ),\
            "Data/MeanFieldGS.dat" u 1:2 w l ls 11 title " MeanField"
}

# plot with same T
do for [jj = 1:3] {

    outname = sprintf("Plots/GroundStates_T%s.pdf", word(TLABELS,jj))  # Note: check that the folder "Plots" exists before running
    
    # save the picture
    set term pdf color enhanced size 5,4
    set output outname
    plot for [idx = 1:9:3] word(FILES,idx-1+jj) u 3:4 w l ls (idx+3)/3 title sprintf( "N=%s; T=%s", word(NLABELS,idx), word(TLABELS,jj) ),\
            "Data/MeanFieldGS.dat" u 1:2 w l ls 11 title " MeanField"
}

set style line  1 dt 2 lw 2 lc "blue"
set style line  2 dt 4 lw 2 lc "web-green"
set style line  3 dt 5 lw 2 lc "red"

# setting table to subtract meanfield results
set table $meanfield
plot "Data/MeanFieldGS.dat" using 2 with table
unset table

set xlabel "{/Symbol l}"; set ylabel "[e_{RSRG}({/Symbol l}) - e_{MF}({/Symbol l})] "


# plot difference with same N and MeanField
do for [jj = 1:9:3] {

    outname = sprintf("Plots/GS_diff_N%s.pdf", word(NLABELS,jj))  # Note: check that the folder "Plots" exists before running
    
    # save the picture
    set term pdf color enhanced size 5,4 dashed
    set output outname       
    plot for [idx = 1:3] word(FILES,idx-1+jj) u 3:($4-$meanfield[$0+1]) w l ls idx title sprintf("N=%s; T=%s", word(NLABELS,jj), word(TLABELS,idx))
}

# plot difference with same T and MeanField
do for [jj = 1:3] {

    outname = sprintf("Plots/GS_diff_T%s.pdf", word(TLABELS,jj))  # Note: check that the folder "Plots" exists before running
    
    # save the picture
    set term pdf color enhanced size 5,4 dashed
    set output outname
    plot for [idx = 1:9:3] word(FILES,idx-1+jj) u 3:($4-$meanfield[$0+1]) w l ls (idx+3)/3 title sprintf("N=%s; T=%s", word(NLABELS,idx), word(TLABELS,jj))
}

