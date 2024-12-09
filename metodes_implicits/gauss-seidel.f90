
! Here we'll create a subroutine (pasitoGS) to calculate the kth step in the Gauss-Seidel method
MODULE gauss_seidel_pas
    IMPLICIT NONE (type, external)

    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)

    INTEGER :: i, j, rows, cols
    REAL(KIND=DP) :: suma1, suma2, old

CONTAINS
SUBROUTINE pasitoGS(mat,vec,b,bigdif)
    REAL(KIND=DP), INTENT(IN) :: mat(:,:) , b(:)
    ! mat is the matrix of the system
    ! b is the vector of the constant terms of the system
    REAL(KIND=DP), INTENT(INOUT) :: vec(:)
    ! vec is the initial vector (step k-1); however, it will be modified and
    ! become the vector step k, which is what the subroutine will return
    REAL(KIND=DP), INTENT(OUT) :: bigdif
    ! We'll use bigdif to keep track of the biggest difference between steps

    bigdif=0
    ! Apply Gauss-Seidel method formula
    DO i=1,rows
        suma1=0
        suma2=0
        DO j=i+1,cols
            suma1=suma1+mat(i,j)*vec(j)
        END DO
        DO j=1,i-1
            suma2=suma2+mat(i,j)*vec(j)
        END DO
        old=vec(i)
        vec(i) = (1/mat(i,i))*( b(i) - suma1 - suma2  )

        ! Store in bigdif the biggest difference of components of vec between steps
        bigdif = MAX(bigdif, ABS(old - vec(i)))

    END DO
END SUBROUTINE pasitoGS
END MODULE gauss_seidel_pas

! Now we'll create a subroutine for the whole method.
MODULE gauss_seidel_tol
    USE gauss_seidel_pas
    IMPLICIT NONE (type, external)

    INTEGER :: numsteps
    REAL(KIND=DP) :: error
    ! We'll use these variables later

CONTAINS
SUBROUTINE gauseideltol(mat,vec,b,tolerance,stepslimit,print)
    ! We'll make the iteration stop when it reaches the tolerance we want (we'll call it tolerance)
    ! However, when the method diverges we want to make sure it stops at some point (and warns us), so we'll
    ! only allow it to take a certain maximum number of steps (we'll call it stepslimit)
    REAL(KIND=DP), INTENT(IN) :: mat(:,:) , b(:) , tolerance
    ! mat and b are the same as before
    INTEGER, INTENT(IN) :: stepslimit
    LOGICAL, INTENT(IN) :: print
    ! We can input it to be true if we want to print convergence of method, false if we don't want to
    REAL(KIND=DP), INTENT(INOUT) :: vec(:) 
    ! vec is the initial condition vector
    ! it will be modified and become the solution (it's the vector the subroutine will return)

    ! Make sure sizes of matrix and vectors are compatible
    rows = SIZE(mat, 1)
    cols = SIZE(mat, 2)
    IF (SIZE(vec) /= rows) THEN
        PRINT *, "Error: Initial condictions vector size must match the number of rows of the matrix."
        STOP
    ELSE IF (SIZE(b) /= rows) THEN
        PRINT *, "Error: Constant vector size must match the number of rows of the matrix."
        STOP
    END IF

    numsteps=0 ! We'll count how many steps are needed and keep track here
    ! Apply Gauss-Seidel as long as error is bigger than tolerance
    DO WHILE (numsteps .LT. stepslimit)
        CALL pasitoGS(mat,vec,b,error)
        ! Starting from the initial vector, we'll modify vector vec applying the method
        ! In the variable 'error' we'll store the value of bigdif for each step and then compare it
        ! to the tolerance in order to know when to stop the iteration.
        numsteps=numsteps+1

        IF (print) THEN
            IF (error .EQ. 0) THEN      ! This is the solution
                WRITE(*,'(A,1X,I0,1X,A)') 'The method converged in ', numsteps , ' iterations. The solution is: '
                DO i=1,rows
                    WRITE(*,*) '$x_{',i,'} = $', vec(i)
                END DO
                EXIT
            ELSE IF (error .LT. tolerance) THEN   ! Stop if we've reached the precision we want
                WRITE(*,*) 'The method converges. The desired precision was reached after ', &
                numsteps, ' steps. The approximate solution is: '
                DO i=1,rows
                    WRITE(*,*) '$x_{',i,'} = $', vec(i)
                END DO
                EXIT
            END IF
        END IF

    END DO

    IF (numsteps+1 .EQ. stepslimit) THEN
        WRITE(*,*) 'The method does not converge.'
    END IF

END SUBROUTINE gauseideltol
END MODULE gauss_seidel_tol
