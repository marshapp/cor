Benvinguts al nostre repositori!

Autors: Núria Castillo Ariño, Eira Jacas García, Marcel López Freixes i Isaac Baldi Garcia

Els programes corresponents al mètode d'Euler explícit es troben a la carpeta "euler_explicit". Els corresponents als mètodes d'Euler implícit i Crank-Nicolson són a la carpeta "metodes_implicits".

Presentem la solució al problema en l'arxiu "solucio.f90", que una vegada compilat i executat hauria d'imprimir a la terminal el resultat del problema (temps màxim que podem aplicar el senyal de 40 V sense violar les condicions que garanteixen l'eficiència del tractament) en segons.

Podem trobar també l'arxiu "errors.gpi", que una vegada carregat mostra dos gràfics amb els errors numèrics de cada mètode: un en el que es comparen tots els mètodes amb tots els mallats proposats ("errors.png") i un altre en el que es comparen tots els mètodes amb un mateix mallat que utilitzem en l'anàlisi de resultats de l'informe pdf ("errors_25.png"). Ara bé, aquest arxiu s'ha d'executar una vegada s'hagin executat els programes relatius al mètode d'Euler explícit, d'Euler implícit (per Gauss-Seidel) i de Crank-Nicolson (per Gauss-Seidel) i conseqüentment generat els arxius amb les dades necessàries següents: "error_explicit.dat", "error_implicit.dat" i "error_crank.dat".

També presentem una animació del probelma a la carpeta "EXTRA_animacio".

Les instruccions referents a l'execució del codi (els READMEs) del mètode d'Euler explícit, dels mètodes implícits (Euler implícit i Crank-Nicolson) i de l'animació són a les seves respectives carpetes.


