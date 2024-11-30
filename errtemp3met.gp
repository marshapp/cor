set terminal png

set key bottom right

set xrange[0:0.02]
set yrange[0:0.004]

set output 'errtemp3met.png'
set xlabel 'z(m)/L'
set ylabel 'Error absolut: ΔT (°C)'
set title 'Error Numèric Euler explícit, Euler implícit i Cranl Nicolson γ=0.25'
plot 'cranc_error025.txt' using 1:2 with lines linecolor rgb "blue" title 'Crank-Nicolson', "resultats_error_v2.txt" u 1:3 with lines linecolor rgb "red" title "Euler explícit" 
#"resultats_error_implicit.txt" u 1:3 with lines linecolor rgb "green" title "Euler implícit"