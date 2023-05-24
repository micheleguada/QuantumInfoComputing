##### Fit of spacing distribution #####

datafile = ARG1
Xmax = ARG2
Ymax = ARG3
length = strlen(datafile)
MatType = substr(datafile, 0, length-4)

## Defining model and fitting
PP(x) = a*(x**alpha)*exp(-b*(x**beta))
fit PP(x) datafile via a, alpha, b, beta

## Set plot options
#set style fill solid
set key right top
set samples 1000
set yrange [0:Ymax]
set xrange [0:Xmax]
set xlabel "Spacing s"; set ylabel "Probability Distribution P(s)"
set grid

## Set labels with parameters
set style textbox opaque noborder
set obj 10 rect from Xmax*0.66,0.9*Ymax to Xmax*0.95,Ymax*0.6
set label 1 "model: P(s) = as^{/Symbol a}exp(-bs^{/Symbol b})" at Xmax*0.68,Ymax*0.85 boxed font "arialbd,14"
set label 2 sprintf("a = %3.4f", a)               at Xmax*0.71,Ymax*0.80 boxed font "arialbd,14"
set label 3 sprintf("{/Symbol a} = %3.4f", alpha) at Xmax*0.71,Ymax*0.75 boxed font "arialbd,14"
set label 4 sprintf("b = %3.4f", b)               at Xmax*0.71,Ymax*0.70 boxed font "arialbd,14"
set label 5 sprintf("{/Symbol b} = %3.4f", beta)  at Xmax*0.71,Ymax*0.65 boxed font "arialbd,14"

## Save the plot
set term pdf color enhanced size 7,5
set output sprintf("%s%s%s","Fit_",MatType,".pdf")
plot datafile w boxes lc "red" title MatType, PP(x) w l lc "blue" lw 1.2 title sprintf("%s%s", MatType, " fit")

## Save Residuals plot
unset yrange
set xlabel "Spacing s"; set ylabel sprintf("%s%s", "Fit Residuals of ",MatType)
set term pdf color enhanced size 7,5
set output sprintf("%s%s%s","Res_",MatType,".pdf")
plot datafile using 1:($2-PP($1)) w lp pt 13 ps 0.7 lc "red" t "Residuals"
