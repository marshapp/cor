PROGRAM crank

    USE gauss_seidel_tol, ONLY: gauseideltol
    
    implicit none (type, external)
    
    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
    
    
    integer, parameter :: N=101 , cv = 3686 , rho = 1081 , V = 40
    real(KIND=DP), parameter :: ta = 0.025, kk = 0.56 , sigma = 0.472 , Tc = 36.5, L = 0.02
    integer :: i , j , m , s , rows, steps, scale , data , error_data , scales(3)
    real(KIND=DP) :: cond , alpha , beta , pext , gamma , dz , dt , position , mat(1:N-2,1:N-2) , bp(1:N-2), tolerance
    real(KIND=DP), allocatable :: ter(:,:), vecm(:), temperaturesta(:,:)
    ! We'll store all the temperature values in ter.
    ! Every row is a certain time and every column a certain position.
    ! Each row is the vector we'll apply the Gauss-Seidel method to in order to find the next row.
    
    

    ! Define numerical constants from the problem.
    pext = sigma*V**2/(2*L**2)
    alpha = kk/(cv*rho)
    beta = pext/(cv*rho)
    cond = Tc*alpha/(beta*L**2)  ! This is Tc normalised


    tolerance = 0.001*alpha/(beta*L**2)  ! Tolerance we want in Gauss-Seidel method.
    ! It's defined by a difference of temperature (in °C) of 0,001.
    steps = 1000 ! ! Maximum number of steps we allow Gauss-Seidel to take.
    
    
    scales=[1,2,4]
    ! We'll use scale to define the relation between dz and dt (see below).
    ! scale is the inverse of gamma.
    ! We keep in scales the values of scale (and therefore gamma) over
    ! which we want to iterate (we'll use them in the report).

    ALLOCATE(vecm(N-2))
    ALLOCATE(temperaturesta(3,N))
    ! In temperatures ta we'll write the temperatures we obtained at ta.
    ! It has one row for every value of scale.
    
    DO s=1,3
        ! Define the discretizations.
        scale=scales(s)
        dz= 1.0_DP / (N - 1)
        dt=(1.0_DP/scale)*dz**2
        gamma = dt/(dz**2)
    
        ! Define the matrix mat, which is the matrix that contains the coefficients of the system of equations.
        DO i = 1,N-2
            mat(i, i) = gamma + 1  ! diag
            IF (i .GE. 2) THEN
                mat(i, i-1) = -gamma/2   ! lower diag
            END IF
            IF (i .LE. N-3) THEN
                mat(i, i+1) = -gamma/2   ! upper diag
            END IF
        END DO
    
        ! Define number of temporal steps we'll take. Depends on dt, and therefore on scale.
        rows = int(ta/dt)+1
        ! According to how we previously defined it, this is what the size of ter should be.
        ALLOCATE(ter(rows,N))
    
        ! Fill temperatures known by CC.
        DO m = 1, rows
            ter(m, 1) = cond
            ter(m, N) = cond
        END DO
    
        ! Define the initial vector which is row 1, ie at t=0, known by IC.
        DO j = 1, N
            ter(1,j) = cond
        END DO
    
    
        ! Find temperatures using Gauss-Seidel method, for each row (ie for each time).
        DO m = 2,rows
            ! Define the components of bp, which is the vector containing the system's constants.
            DO i = 2,N-3
                bp(i) = dt + (gamma/2)*ter(m-1,i+2) + (1-gamma)*ter(m-1,i+1) + (gamma/2)*ter(m-1,i)
            END DO
            ! Since the temperatures at j=1 and j=N are not variables, the equation used to
            ! find the temperature at j=2 and j=N-2 are modified, so their bps are different from the rest.
            bp(1) = dt + (gamma/2)*ter(m-1,3) + (1-gamma)*ter(m-1,2) + (gamma/2)*ter(m-1,1) + (gamma/2)*cond
            bp(N-2) = dt + (gamma/2)*ter(m-1,N) + (1-gamma)*ter(m-1,N-1) + (gamma/2)*ter(m-1,N-2) + (gamma/2)*cond
    
            vecm = ter(m-1,2:N-1) ! We'll use the vector from the previous step as the
            ! initial vector for the iteration, since it's probably close to the one we want to find.

            ! We don't want to print that the method converges every time so we input .false.

            CALL gauseideltol(mat,vecm,bp,tolerance,steps,.false.)  

            ! Define the next row as the result of Gauss-Seidel
            ter(m,2:N-1) = vecm
    
        END DO
    
        temperaturesta(s,:) = ( ter(rows,:)*(beta*(L**2)) )/alpha
        ! We'll store here the values of the temperature (in °C) at ta (which corresponds to last row of ter) for each gamma.
    
        DEALLOCATE(ter)
    
    END DO
    
    
    ! Create a .txt file in which to put the data we need for the graph.
    OPEN(NEWUNIT=data, FILE='crank.dat', STATUS='unknown', ACTION='WRITE')
    
    ! We'll write data for ta for each position.
    ! 1st column: positions (z) ; last column: analytical solution; middle columns: numerical solutions
    ! We also undo the normalisation
    position = 0
    
    DO i = 1,N
        WRITE(data,*) position , temperaturesta(1,i) , temperaturesta(2,i) , temperaturesta(3,i) , f(200,position/L)
        position = position+dz*L
    END DO
    
    
    !----------------- Error ----------------

    ! Create a .txt file in which to put the errors for each position, calculated as the difference between the analytical and the numerical solutions.
    open(NEWUNIT=error_data, file='error_crank.dat', status='unknown', action='write')
        do j = 0, N-1
            write(error_data, *) (j*dz)*L , abs( temperaturesta(1,j+1)-f(200, j*dz) ) , &
            abs( temperaturesta(2,j+1)-f(200, j*dz) ), abs( temperaturesta(3,j+1)-f(200, j*dz) )
        end do
    close(error_data)
    
    
    DEALLOCATE(vecm)
    DEALLOCATE(temperaturesta)


    ! Define the function corresponding to the analytic solution at ta
    ! Its inputs are the normalised position and the number of terms of the summation in its expression
    
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
        f=(f_norm*beta*(L**2))/alpha ! Denormalise the result
    
    end function
    
    END PROGRAM
