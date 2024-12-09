
Aquesta carpeta conté tots els programes que hem usat per generar les dues animacions gif que presentem al informe pdf i ambdues animacions ja generades.
Les animacions són als fitxers: animacio1.gif i animacio2.gif .

El fitxer .f90 conté el codi de fortran per generar les dades que necessitem per fer l'animació. 
En ell es soluciona l'equació diferencial del nostre problema per uns quants punts temporals representatius i es generen els fitxers .txt corresponents a cada temps amb els valors de les temperatures a cada punt de l'espai. En aquest cas usem el mètode de Crank Nicolson ja que no necessitem dades amb molta precissió.

El fitxer .gp conté el codi de gnuplot per passar les dades dels fitxers .txt a imatges .png que seràn els frames del gif. 
Si es vol obtenir la primera animació s'ha de desactivar l'opció $set palette maxcolors 7 i així obtindrem un degradat de colors. 
La segona animació es genera activant aquesta opció.

Per últim per generar el gif, usem ImageMagick. A la terminal escrivim el següent codi:
$magick convert -delay 10 -loop 0 *.png animacio.gif
Amb  el qual compilem les imatges en un sol arxiu gif.

