Primerament s'ha de compilar i executar el programa "jacobi_implicit.f90". Un cop fet això, s'haurien de generar automàticament dos arxius de dades que seran necessaris per fer els gràfics: "jacobi_implicit.dat" i "error_jacobi_implicit.dat".

Tanmateix, abans de generar els gràfics cal haver compilat primer el fitxer "implicit.f90" (seguint les instruccions del fitxer README.md corresponent), ja que s'utilitzen les dades obtingudes d'aquest ("implicit.dat" i "error_implicit.dat") en els gràfics d'aquesta secció.

Amb això, si es carrega "jacobi_resultats.gpi", haurien de generar-se automàticament dues imatges: "jacobi_implicit_comparacio.png" i "error_jacobi_implicit_comparacio.png". En aquests dos gràfics es comparen els resultats obtinguts pel mètode de Gauss-Seidel i el de Jacobi, utilitzats per implementar el mètode d'Euler implícit amb gamma=1.
