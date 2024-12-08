És important seguir aquests passos per una execució correcta del nostre codi:

Pel que fa al mètode explicit, primer s'ha de compilar a la terminal "explicit.f90" i un cop fet això s'haurien de crear dos arxius de dades necessaris per fer els gràfics: "explicit.dat" i "error_explicit.dat". Amb això, si es carrega "explicit.gpi", haurien de crear-se dues imatges: "explicit_51.png" i "explicit_49i25.png", amb les solucions graficades pels diferents mallats.

Pel que fa als mètodes implicits, a la corresponent carpeta s'hi trobaran els arxius necessàris, sent el mètode de Gauss-Sidel el nostre mètode principal de resolució! Presentem la resolució per jacobi com un EXTRA a mode de comparació per completitud. El cas és que per la solució numèrica pel mètode d'Euler implicit i Crank-Nicolson utilitzem un subprograma que conté un mòdul amb el mètode de gauss-seidel per la resolució de sistemes d'equacions, per tant, cal posar atenció a l'execució d'aquesta part del codi. Cal seguir estrictament el següent ordre:

-Primer hem de compilar el subprograma que conté el mòdul referent a gauss-seidel, "gauss-seidel.f90", en un arxiu objecte escribint a la terminal "gfortran gauss-seidel.f90 -o gauss-seidel.o -c".
-Ara, una vegada fet això, per compilar els programes corresponents al mètode d'Euler implicit i de Crank-Nicolson, s'ha d'escriure a la terminal el següent: Primer, "gfortran -c implicit.f90 -o implicit.o" i després "gfortran implicit.o -o implicit" i el mateix amb crank.f90.

Una vegada s'executin els programes amb "./implicit" o "./crank" es crearran dos arxius .dat: "implicit.dat" i "error_implicit.dat" o "crank.dat" i "error_crank.dat" que permetran la creació de les imatges que contenen els gràfics de les solucions corresponents carregant "implicit.gpi" o "crank.gpi".


