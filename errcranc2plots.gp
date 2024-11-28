set terminal png

set key bottom right

set style line 1 lc rgb "red" lw 2 lt 1 pt 7 ps 1
set style line 2 lc rgb "blue" lw 2 lt 1 pt 7 ps 1

set xrange[0:0.02]
set yrange[0:0.0035]


set output 'errcranc2plots.png'
set xlabel 'z(m)/L'
set ylabel 'error absolut'
set title 'Error Numèric Crank Nicolson'
plot 'cranc_error1.txt' using 1:2 with lines title 'gamma = 1.0', 'cranc_error05.txt' using 1:2 with lines title 'gamma = 0.5'