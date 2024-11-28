set terminal png
set output 'errcranc.png'
set xlabel 'z(m)/L'
set ylabel 'error absolut'
set title 'Error Numèric Absolut Cranc Nicolson'
plot 'cranc_error.txt' using 2:5 with lines title 'Crank 0.25'
# 5 es l'error absolut a l'eix y, 2 per a z/L i 1 per a z en l'eix x.