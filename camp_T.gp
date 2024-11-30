set terminal pngcairo size 800,600 enhanced font "Arial,14" 
set output 'camp.png'

set multiplot
set autoscale fix
set colorbox
set tics
#set palette defined ( 0 "black", 1 "red", 2 "yellow" )
set palette defined ( \
    36.5 "black", \
    40 "red", \
    45 "yellow", \
    51 "yellow",\
    51 "blue",\
    55 "blue" )
set cblab 'Temperatura (ºC)'
set xlabel 'z(m)'
set ylabel 't(s)'
#set parametric
set title 'Camp de temperatura al llarg del temps'
#set arrow from 0.007525,graph 0 to 0.007525,graph 1  nohead lw 5 lc rgb "black"
plot 'camp.txt' with image notitle

set rmargin at screen 1089.0
set xrange [0:0.02]
set yrange [0:0.025]
set autoscale fix
set colorbox
set palette defined ( 0 "black", 1 "red", 2 "yellow" )
set cblab 'Temperatura (ºC)'
unset border
unset key
unset xlabel
unset ylabel
unset title
unset tics
plot 'linia1.txt' with lines lw 2 lc 'green', 'linia2.txt' with lines lw 2 lc 'green'



