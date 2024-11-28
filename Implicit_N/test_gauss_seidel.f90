PROGRAM TestGaussSeidel
    !USE gauss_seidel_ss, ONLY: gauseidelss
    USE gauss_seidel_tol, ONLY: gauseideltol

    IMPLICIT NONE (type, external)

    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)

    INTEGER, PARAMETER :: N = 3
    REAL(KIND=DP) :: mat(N, N), vec(N), b(N)
    REAL(KIND=DP) :: tolerance
    INTEGER :: maxsteps

    ! Initialize matrix, vector, tolerance and maxsteps
    b= [1,-1,1] !escriu el vector del terme independent, b    
    mat= reshape([11,-1,1, &
                 -1,6,1, &
                 1,1,7], [n, n])!definim la matriu A

    vec= [0.0, 0.0, 0.0] ! Initial guess
    tolerance = 1.0e-7
    maxsteps = 50

    ! Call Gauss-Seidel
    
    !WRITE(*,*) 'GS ss'
    !CALL gauseidelss(mat, vec, b, tolerance, maxsteps, .true.)
    WRITE(*,*) 'GS tol'
    CALL gauseideltol(mat, vec, b, tolerance,1000,.true.)

END PROGRAM TestGaussSeidel