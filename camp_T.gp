set terminal png
set output 'camp.png'
set autoscale fix
set xtics 1
set ytics 1
set palette defined ( 0 "black", 1 "red", 2 "yellow" )
set cblab 'Temperatura (ºC)'
set xlabel 'z(m)'
set ylabel 't(s)'
set title 'Camp de temperatura al llarg del temps'

plot 'camp.txt' with image notitle