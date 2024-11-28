set terminal png
set output 'errcranc.png'
set xlabel 'z(m)/L'
set ylabel 'error absolut'
set title 'Error Numèric Crank Nicolson'
plot 'cranc_error.txt' using 1:2 with lines title 'Crank 1'
# 5 es l'error absolut a l'eix y, 2 per a z/L i 1 per a z en l'eix x