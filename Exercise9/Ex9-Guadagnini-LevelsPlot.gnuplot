# gnuplot script to plot the first kk eigenvalues in function of Lambda
#       for different N
# To run: gnuplot Ex9-Guadagnini-LevelsPlot.gnuplot 

# palette: paired
# line styles
set style line 1 lt 1 lw 2 pt 7 ps 0.3 lc rgb '#A6CEE3' # light blue
set style line 2 lt 1 lw 2 pt 7 ps 0.3 lc rgb '#1F78B4' # dark blue
set style line 3 lt 1 lw 2 pt 7 ps 0.3 lc rgb '#B2DF8A' # light green
set style line 4 lt 1 lw 2 pt 7 ps 0.3 lc rgb '#33A02C' # dark green
set style line 5 lt 1 lw 2 pt 7 ps 0.3 lc rgb '#FB9A99' # light red
set style line 6 lt 1 lw 2 pt 7 ps 0.3 lc rgb '#E31A1C' # dark red
set style line 7 lt 1 lw 2 pt 7 ps 0.3 lc rgb '#FDBF6F' # light orange
set style line 8 lt 1 lw 2 pt 7 ps 0.3 lc rgb '#FF7F00' # dark orange

# graphical options
set xlabel "{/Symbol l}"; set ylabel "E_n({/Symbol l}) /(N-1)"
set key bottom left
set grid
# set format y "%.3g"
set pointsize 0.5

k_in = 8 #5
# Loop over different N
do for [NN = 2:12] {  
    # build the names
    in_name = sprintf("Data/EigVals_N%.2d.dat",  NN)
    outname = sprintf("Plots/EigVals_N%.2d.pdf", NN) # Note: check that the folder "Plots" exists before running
    
    if (k_in > 2**NN) { kk = 2**NN } else { kk = k_in }
    
    set style textbox opaque
    set label 1 sprintf("N=%d",NN) at graph 0.5,0.05 center boxed font "arialbd,14"
    
    # save the picture
    set term pdf color enhanced size 5,4
    set output outname
    plot for [i=1:kk] in_name u 1:(column(1+i)/(NN-1)) w l ls i title "E_".(i-1)."({/Symbol l})"
}

# Loop over eigenvalues for different N
do for [kk = 1:8] { 
    # build the names
    outname = sprintf("Plots/EigVals_k%.2d.pdf", (kk-1)) # Note: check that the folder "Plots" exists before running
 
    set key default
    set ylabel "E_".(kk-1)."({/Symbol l}) / (N-1)"
    unset label
    
    set palette defined (0 'red',  1 'blue')
    set cbrange [3:12]
    set cblabel "N"
    
    # save the picture
    set term pdf color enhanced size 5,4
    set output outname
    plot for [NN=3:12] sprintf("Data/EigVals_N%.2d.dat",NN) u 1:(column(1+kk)/(NN-1)) w l lt palette frac (NN-3)/10. lw 1.5 title ""
}

