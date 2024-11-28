PROGRAM IVEGOTTHIS

USE gauss_seidel_tol, ONLY: gauseideltol

implicit none (type, external)

INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)

integer, parameter :: N=101 , cv = 3686 , rho = 1081 , V = 40
real(KIND=DP), parameter :: ta = 0.025, kk = 0.56 , sigma = 0.472 , Tc = 36.5, L = 0.02
integer :: i , j , m , rows, steps, scale , data
real(KIND=DP) :: cond , alpha , beta , pext , gamma , dz , dt , position , mat(1:N-2,1:N-2) , bp(1:N-2), tolerance
real(KIND=DP), allocatable :: ter(:,:), vecm(:), temperaturesta(:)
! We'll store all the temperature values in ter.
! Every row is a certain time and every column a certain position.
! Its rows are the vectors we'll input in gauss-seidel.

scale=1 ! Scale defines relation between dz and dt (see below).
tolerance = 10**(-3) ! Tolerance we want in G-S method.
steps = 5000 ! Max number of steps we allow G-S to take.

dz= 1.0_DP / (N - 1)  ! Use 1.0_DP to force floating-point division.
dt=(1.0_DP/scale)*dz**2

pext = sigma*V**2/(2*L**2)
alpha = kk/(cv*rho)
beta = pext/(cv*rho)
cond = Tc*alpha/(beta*L**2)
gamma = dt/(dz**2)

! Define the matrix mat.
DO i = 1,N-2
    mat(i, i) = 2 * gamma + 1  ! diag
    IF (i .GE. 2) THEN
        mat(i, i-1) = -gamma   ! lower diag
    END IF
    IF (i .LE. N-3) THEN
        mat(i, i+1) = -gamma   ! upper diag
    END IF
END DO


! Define size of ter as well as allocatable vectors.
rows = CEILING(ta/dt)
ALLOCATE(ter(rows,N))
ALLOCATE(vecm(N-2))
ALLOCATE(temperaturesta(N))

! Known by CC.
DO m = 1, rows
    ter(m, 1) = cond
    ter(m, N) = cond
END DO

! Define the initial vector which is row 1, ie at t=0, known by IC.
DO j = 1, N
    ter(1,j) = cond
END DO



! CHECKS
write(*,*) 'dz=',dz,'dt=',dt,'gamma=',gamma

WRITE(*,*) 'mat'
DO i=1,6
    WRITE(*,*) mat(i,1:6)
END DO

WRITE(*,*) 'ter'
DO i=1,6
    WRITE(*,*) ter(i,1:3),ter(i,N-1:N)
END DO




! Find temperatures using method.
DO m = 2,rows
    ! Define the components of bp.
    DO i = 2,N-3
        bp(i) = dt + ter(m-1,i)
    END DO
    bp(1) = dt + cond*(gamma) + ter(m-1,2)
    bp(N-2) = dt + cond*(gamma) + ter(m-1,N-1)

    vecm = ter(m-1,2:N-1) ! We'll use the vector from the previous step as the
    ! starting point, since it's probably close to the one we want to find.

    CALL gauseideltol(mat,vecm,bp,tolerance,steps,.false.)

    ter(m,2:N-1) = vecm

END DO

! CHECK
WRITE(*,*) 'ter resultant'
DO i = 1,7
    WRITE(*,*) ter(i,1:3)
END DO
WRITE(*,*) 'ultima fila ter = ',  ter(N,1:4)

! Create a .txt file in which to put the data we need for the graph.
OPEN(NEWUNIT=data, FILE='data_implicit_gs.dat', STATUS='unknown', ACTION='WRITE')

! We'll write data for t=ta (which corresponds to last row of ter)
! 1st column: positions (z) ; 2nd column: temperatures
! We also undo the normalisation
position = 0
temperaturesta = ( ter(rows,:)*(beta*(L**2)) )/alpha

write(*,*) 'ta(1:6): ', temperaturesta(1:6)

DO i = 1,N
    WRITE(data,*) position , temperaturesta(i)
    position = position+dz*L
END DO

CLOSE(data)



DEALLOCATE(ter)
DEALLOCATE(vecm)
DEALLOCATE(temperaturesta)
END PROGRAM