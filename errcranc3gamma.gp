set terminal png

set key bottom right

set xrange[0:0.02]
set yrange[0:0.004]

set output 'errcranc3gamma.png'
set xlabel 'z(m)/L'
set ylabel 'Error absolut: ΔT (°C)'
set title 'Error Numèric Crank Nicolson'
plot 'cranc_error1.txt' using 1:2 with lines linecolor rgb "blue" title 'gamma = 1.0', 'cranc_error05.txt' using 1:2 with lines linecolor rgb "red" title 'gamma = 0.5', 'cranc_error025.txt' using 1:2 with lines linecolor rgb "green" title 'gamma = 0.25'