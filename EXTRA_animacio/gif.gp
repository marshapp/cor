
set terminal pngcairo size 800,600 enhanced font "Arial,14"


# Creem els frames del fig des dels fitxers .txt amb un bucle
#Usem sprintf per generar el nom de cada fitxer de sortida i el títol amb el temps corresponent a cada frame
#Usem maxcolors perquè ens quedi un camp de temperatures amb colors plans, si volem degradats s'ha de desactivar aquesta opció
#usem multiplot per combinar el camp de temperatures amb les marques limitadores del teixit biològic dolent
#usem set palette defined per definir els colors de les diferents temperatures
do for [frame=1:101] {
    reset
    set output sprintf("frame_%03d.png", frame)
    set multiplot
    set autoscale fix
    set border
    set title
    set tics
    set colorbox
    set cblab 'Temperatura (ºC)'
    set xlabel 'x(m)'
    set ylabel 'y(m)'
    set title sprintf("Camp de temperatura al temps: %d s", frame*0.63)

    set palette maxcolors 7
    set palette defined ( \
        36.5 "black", \
        40 "brown", \
        45 "red", \
        47 "orange",\
        50 "dark-yellow",\
        52 "yellow",\
        60 "white" )
    set cbrange[36.5:60]
    
    
    plot sprintf("d_%d.txt", frame) with image notitle
    
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
    plot 'linia1.txt' with lines lw 2 lc 'black', 'linia2.txt' with lines lw 2 lc 'black'
    unset multiplot
    unset xrange 
    unset yrange
    
}
