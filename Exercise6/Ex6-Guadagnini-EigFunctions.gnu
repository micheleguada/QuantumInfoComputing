##### Plot of first 4 eigenfunctions #####

## Set plot options
set key right top
#set yrange [0:Ymax]
set xrange [-8:8]
set xlabel "x"; set ylabel "{/Symbol y}_i(x)"
set grid

# palette
set palette maxcolors 8
set palette defined ( 0 '#1B9E77',\
    	    	      1 '#D95F02',\
                      2 '#7570B3',\
                      3 '#E7298A',\
                      4 '#66A61E',\
                      5 '#E6AB02',\
                      6 '#A6761D',\
                      7 '#666666' )

## Save the plot
set term pdf color enhanced size 6,4
set output "Eigenfunctions.pdf"
plot "Eigenfunctions.dat" u 1:2 w l lw 1.5 title "{/Symbol y}_1(x),  E_1=1",\
     "Eigenfunctions.dat" u 1:3 w l lw 1.5 title "{/Symbol y}_2(x),  E_2=3",\
     "Eigenfunctions.dat" u 1:4 w l lw 1.5 title "{/Symbol y}_3(x),  E_3=5",\
     "Eigenfunctions.dat" u 1:5 w l lw 1.5 title "{/Symbol y}_4(x),  E_4=7"
     
