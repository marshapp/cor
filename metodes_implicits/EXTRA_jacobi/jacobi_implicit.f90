PROGRAM jacobi
    
    implicit none (type, external)
    
    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
    
    integer, parameter :: N=101 , cv = 3686 , rho = 1081 , V = 40
    real(KIND=DP), parameter :: ta = 0.025, kk = 0.56 , sigma = 0.472 , Tc = 36.5, L = 0.02
    integer :: i , j , m , rows, steps, scale , data , error_data
    real(KIND=DP) :: cond , alpha , beta , pext , gamma , dz , dt , position , mat(1:N-2,1:N-2) , bp(1:N-2), tolerance
    real(KIND=DP), allocatable :: ter(:,:), vecm(:), temperaturesta(:)
    ! We'll store all the temperature values in ter.
    ! Every row is a certain time and every column a certain position.
    ! Its rows are the vectors we'll input in gauss-seidel.


    !Define the variables to run jacobi method
    !in x vector we'll store the solution of each iteration
    !x_1 contains the previous solution
    real(KIND=DP) ::  x(N-2), x_1(N-2)
    integer :: iter
    real(KIND=DP) :: sum, error


    ! Define numerical constants from the problem.
    pext = sigma*V**2/(2*L**2)
    alpha = kk/(cv*rho)
    beta = pext/(cv*rho)
    cond = Tc*alpha/(beta*L**2)   ! This is Tc normalised.


    tolerance = 0.001*alpha/(beta*L**2)  ! Tolerance we want in Jacobi method.
    ! It's defined by a difference of temperature (in °C) of 0,001.
    steps = 1000 ! Maximum number of steps we allow'll Gauss-Seidel to take.
    
    scale=1
    ! We'll use scale to define the relation between dz and dt (see below).
    ! scale is the inverse of gamma.
    ! We only do it for one value of scale because we'll only compare this particular case.

    ! Define the discretizations.
    dz= (1.0_DP) / (N - 1)
    dt=(1.0_DP/scale)*(dz**2)
    
    
    ! Define the matrix mat,  which is the matrix that contains the coefficients of the system of equations.
    DO i = 1,N-2
        mat(i, i) = 2 * gamma + 1  ! diag
        IF (i .GE. 2) THEN
            mat(i, i-1) = -gamma   ! lower diag
        END IF
        IF (i .LE. N-3) THEN
            mat(i, i+1) = -gamma   ! upper diag
        END IF
    END DO

    ALLOCATE(vecm(N-2))
    ALLOCATE(temperaturesta(N))
    ! In temperatures ta we'll write the temperatures we obtained at ta.

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
    
    ! Find temperatures using Jacobi method, for each row (ie for each time).
    DO m = 2,rows
        ! Define the components of bp, which is the vector containing the system's constants.
        ! Element i of bp corresponds to element i+1 of ter.
        DO i = 2,N-3
            bp(i) = dt + ter(m-1,i+1)
        END DO
        ! Since the temperatures at j=1 and j=N are not variables, the equation used to
        ! find the temperature at j=2 and j=N-2 are modified, so their bps are different from the rest.
        bp(1) = dt + cond*(gamma) + ter(m-1,2)
        bp(N-2) = dt + cond*(gamma) + ter(m-1,N-1)
    
        

        x= ter(m-1,2:N-1)   ! We'll use the vector from the previous step as the
        ! initial vector for the iteration, since it's probably close to the one we want to find.

        ! Iterate through the Jacobi method.
        do iter = 1, steps
            ! The solution of the temperature vector in the previous step is used to find the next one.
            x_1 = x 
            ! Now we implement Jacobi formula for each component of the temperatures vector.
            do i = 1, N-2 
                sum = 0
                do j = 1, N-2
                    if (j/=i) then
                        sum = sum + mat(i,j)*x_1(j)
                    end if
                end do 
                ! Here we find the new value for the temperature vector.
                x(i)=(bp(i)-sum)/mat(i,i) 
            end do
            ! Compute the difference betwen contiguous solutions (ie error).
            error = maxval(abs(x - x_1))

            ! When error is below tolerance iteration stops.
            if (error < tolerance) then
                exit
            end if
            if (iter==steps) then
                stop 
            end if
            
        end do
    
        ! Store the solution in the ter matrix.
        ter(m,2:N-1) = x
    
    END DO


    
    ! Create a .txt file in which to put the data we need for the graph.
    OPEN(NEWUNIT=data, FILE='jacobi_implicit.dat', STATUS='unknown', ACTION='WRITE')
    
    ! We'll write data for ta for each position.
    ! 1st column: positions (z) ; 2nd column: numerical solutions
    ! We also undo the normalisation
    position = 0
    temperaturesta = (ter(rows,:)*(beta*(L**2)))/alpha
    
    DO i = 1,N
        WRITE(data,*) position , temperaturesta(i)
        position = position+dz*L
    END DO
    
    CLOSE(data)

    DEALLOCATE(vecm)
    DEALLOCATE(ter)

    !----------------- Error ----------------

! Create a .txt file in which to put the error for each position, calculated as the difference between the analytical and the numerical solution.
open(NEWUNIT=error_data, file='error_jacobi_implicit.dat', status='unknown', action='write')
do j = 0, N-1
    write(error_data, *) (j*dz)*L , abs( temperaturesta(j+1)-f(300, j*dz) )
end do
close(error_data)


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
f=(f_norm*beta*(L**2))/alpha   ! Denormalise the result

end function
    
    

END PROGRAM
