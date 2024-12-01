PROGRAM IVEGOTTHIS

USE gauss_seidel_tol, ONLY: gauseideltol

implicit none (type, external)

INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)


integer, parameter :: N=101 , cv = 3686 , rho = 1081 , V = 40
real(KIND=DP), parameter :: ta = 0.025, kk = 0.56 , sigma = 0.472 , Tc = 36.5, L = 0.02
integer :: i , j , u , m , Nf , rows, steps, scale , data , error_data
real(KIND=DP) :: cond , alpha , beta , pext , gamma , dz , dt , position , mat(1:N-2,1:N-2) , bp(1:N-2), tolerance
real(KIND=DP), allocatable :: ter(:,:), vecm(:), temperaturesta(:,:)
! We'll store all the temperature values in ter.
! Every row is a certain time and every column a certain position.
! Its rows are the vectors we'll input in gauss-seidel.



tolerance = (10.0_DP)**(-3) ! Tolerance we want in G-S method.
steps = 5000 ! Max number of steps we allow G-S to take.

pext = sigma*V**2/(2*L**2)
alpha = kk/(cv*rho)
beta = pext/(cv*rho)
cond = Tc*alpha/(beta*L**2)

ALLOCATE(vecm(N-2))
ALLOCATE(temperaturesta(2,N))


DO scale = 1,2 ! Scale defines relation between dz and dt (see below).
    dz= 1.0_DP / (N - 1)  ! Use 1.0_DP to force floating-point division.
    dt=(1.0_DP/scale)*dz**2
    gamma = dt/(dz**2)

    ! Define the matrix mat, which contains the coefficients of the system of equations for each t.
    DO i = 1,N-2
        mat(i, i) = 2 * gamma + 1  ! diag
        IF (i .GE. 2) THEN
            mat(i, i-1) = -gamma   ! lower diag
        END IF
        IF (i .LE. N-3) THEN
            mat(i, i+1) = -gamma   ! upper diag
        END IF
    END DO

    ! Define size of ter.
    rows = CEILING(ta/dt)
    ALLOCATE(ter(rows,N))

    ! Known by CC.
    DO m = 1, rows
        ter(m, 1) = cond
        ter(m, N) = cond
    END DO

    ! Define the initial vector which is row 1, ie at t=0, known by IC.
    DO j = 1, N
        ter(1,j) = cond
    END DO


    ! Find temperatures using Gauss-Seidel method.
    DO m = 2,rows
        ! Define the components of bp.
        DO i = 2,N-3
            bp(i) = dt + ter(m-1,i+1)
        END DO
        bp(1) = dt + cond*(gamma) + ter(m-1,2)
        bp(N-2) = dt + cond*(gamma) + ter(m-1,N-1)

        vecm = ter(m-1,2:N-1) ! We'll use the vector from the previous step as the
        ! starting point, since it's probably close to the one we want to find.

        CALL gauseideltol(mat,vecm,bp,tolerance,steps,.false.)

        ter(m,2:N-1) = vecm

    END DO

    temperaturesta(scale,:) = ( ter(rows,:)*(beta*(L**2)) )/alpha
    ! We'll store here the values of the temperature at t=ta (in °C) for each gamma.

    DEALLOCATE(ter)

END DO


! Create a .txt file in which to put the data we need for the graph.
OPEN(NEWUNIT=data, FILE='data_implicit_gs.dat', STATUS='unknown', ACTION='WRITE')

! We'll write data for t=ta (which corresponds to last row of ter)
! 1st column: positions (z) ; 2nd column: temperatures
! We also undo the normalisation
position = 0

DO i = 1,N
    WRITE(data,*) position , temperaturesta(1,i) , temperaturesta(2,i) , f(300,position/L)
    position = position+dz*L
END DO


!----------------- Error ----------------

open(NEWUNIT=error_data, file='resultats_error_implicit.dat', status='unknown', action='write')
    ! We write in error_data the error for each position (en m)
    do j = 0, N-1
        write(error_data, *) (j*dz)*L , abs( temperaturesta(1,j+1)-f(200, j*dz) ) , abs( temperaturesta(2,j+1)-f(300, j*dz) )
    end do
close(error_data)


DEALLOCATE(vecm)
DEALLOCATE(temperaturesta)

contains

real(KIND=DP) function f(Nf, x_norm)
    implicit none(type, external)
    integer, intent(in) :: Nf
    real(KIND=DP), intent(in) :: x_norm
    real(KIND=DP), parameter :: t_norm=0.025, pi=2*acos(0.0)
    real(KIND=DP) :: sumatori, f_norm
    integer :: k

    sumatori=0 
    do k = 1, Nf
        sumatori = sumatori + ((1-exp(-((2*k-1)**2)*(pi**2)*t_norm))/((2*k-1)**3))*sin((2*k-1)*pi*x_norm)
    end do

    f_norm=cond+(4/(pi**3))*sumatori
    f=(f_norm*beta*(L**2))/alpha

end function

END PROGRAM