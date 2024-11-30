set terminal png

set key bottom right

set xrange[0:0.02]
set yrange[0:0.0035]

set output 'errcranc3met.png'
set xlabel 'z(m)/L'
set ylabel 'error absolut (ºC)'
set title 'Error Numèric Euler explícit, Euler implícit i Cranl Nicolson γ=0,49'
plot 'cranc_error05.txt' using 1:2 with lines title 'Crank-Nicolson', "resultats_error_v2.txt" u 1:2 w l ls 1 t "Euler explícit", "resultats_error_implicit.txt" u 1:3 w l ls 1 t "Euler implícit"