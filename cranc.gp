set terminal png
set output 'cranc.png'
set xlabel 'z(m)'
set ylabel 'T(ºC)'
set title 'Cranc Nicolson'
plot 'cranc_data.txt' using 1:2 with lines title 'Cranc Nicolson'
