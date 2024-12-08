Les instruccions referents a l'execució del codi (els READMEs) del mètode d'Euler explicit i dels mètodes implicits (Euler implicit i Crank-Nicolson) són a les seves respectives carpetes.

L'estructura dels arxius entregats és la següent:

El programa corresponent al mètode d'Euler explicit junt amb el respectiu programa per graficar les diferents solucions es troben a la carpeta "euler_explicit".

Els programes corresponents als mètodes d'Euler implicit i Crank-Nicolson junt amb els seus respectius programes per graficar les solucions són a la carpeta "metodes_implicits". Cal destacar el fet que per fer aquests mètodes hem utilitzat un subprograma amb un mòdul que implementa la solució de sistemes d'equacions per Gauss-Seidel. Ara bé, dins aquesta carpeta, existeix una altra anomenada "EXTRA: jacobi", en aquesta s'hi troben també dos programes que implementen els mètodes d'Euler implicit i jacobi però utilitzant el mètode de jacobi per resoldre els sistemes d'equacions. Ha de quedar clar que la nostra intenció amb aquesta part és la de comparar els resultats obtinguts per Gauss-seidel i per Jacobi. El que pretenem que sigui avaluat principalent és la manera amb la que hem implementat Euler implicit i Crank-Nicolson utilitzant Gauss-Seidel! La part en la que implementem jacobi és extra i per completitud.

A més, presentem també la solució al problema en l'arxiu "solució.f90" que una vegada compilat i executat hauria d'imprimir a la terminal el resultat del problema (temps màxim que podem aplicar el senyal de 40 V sense violar les condicions que garanteixen l'eficiència del tractament) en segons.

Podem trobar també l'arxiu "errors.gpi" que una vegada carregat mostra dos gràfics amb els errors numèrics de cada mètode, un en el que es comparen tots els mètodes amb tots els mallats proposats i un altre en el que es comparen tots els mètodes amb un mateix mallat que utilitzem en l'anàlisi de resultats de l'informe pdf. Ara bé, aquest arxiu s'ha d'executar una vegada s'hagin executat els programes relatius al mètode d'Euler explicit, d'Euler implicit (per Gauss-Sidel) i de Crank-Nicolson (per Gauss-Seidel) i conseqüentment generat els arxius amb les dades necessàries següents: "error_explicit.dat", "error_implicit.dat" i "error_crank.dat".

També presentem un programa per realitzar un gif a la carpeta "EXTRA: gif".


