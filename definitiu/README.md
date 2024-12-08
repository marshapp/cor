L'estructura dels arxius entregats és la següent:

Els programes corresponents al mètode d'Euler explicit es troben a la carpeta "euler_explicit".

Els programes corresponents als mètodes d'Euler implicit i Crank-Nicolson són a la carpeta "metodes_implicits". Cal destacar el fet que per fer aquests mètodes hem utilitzat un subprograma amb un mòdul que implementa la solució de sistemes d'equacions per Gauss-Seidel. 
Dins de la carpeta "metodes_implicits", existeix una altra anomenada "EXTRA: jacobi", en aquesta s'hi troben també programes referents a la implementació del mètode d'Euler implicit però utilitzant el mètode de jacobi per resoldre els sistemes d'equacions. Ha de quedar clar que la nostra intenció és que principalment se'ns avalui la implementació del mètode d'Euler implicit utilitzant Gauss-Seidel i que la part referent a jacobi és un EXTRA que afegim per completitud i en certa manera per justificar l'ús preferent de Gauss-Seidel davant de Jacobi.

Presentem la solució al problema en l'arxiu "solució.f90" que una vegada compilat i executat hauria d'imprimir a la terminal el resultat del problema (temps màxim que podem aplicar el senyal de 40 V sense violar les condicions que garanteixen l'eficiència del tractament) en segons.

Podem trobar també l'arxiu "errors.gpi" que una vegada carregat mostra dos gràfics amb els errors numèrics de cada mètode, un en el que es comparen tots els mètodes amb tots els mallats proposats i un altre en el que es comparen tots els mètodes amb un mateix mallat que utilitzem en l'anàlisi de resultats de l'informe pdf. Ara bé, aquest arxiu s'ha d'executar una vegada s'hagin executat els programes relatius al mètode d'Euler explicit, d'Euler implicit (per Gauss-Sidel) i de Crank-Nicolson i conseqüentment generat els arxius amb les dades necessàries següents: "error_explicit.dat", "error_implicit.dat" i "error_crank.dat".

També presentem una animació del probelma a la carpeta "EXTRA: gif".

A més, tenim aquí també l'informe a entregar com "informe.pdf.

Les instruccions referents a l'execució del codi (els READMEs) del mètode d'Euler explicit, dels mètodes implicits (Euler implicit i Crank-Nicolson) i de l'animació són a les seves respectives carpetes.


