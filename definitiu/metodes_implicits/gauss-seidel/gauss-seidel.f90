! kth step in Gauss-Seidel method
MODULE gauss_seidel_pas
    IMPLICIT NONE (type, external)

    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)

    INTEGER :: i, j, rows, cols
    REAL(KIND=DP) :: suma1, suma2, old

CONTAINS
SUBROUTINE pasitoGS(mat,vec,b,bigdif)
    REAL(KIND=DP), INTENT(IN) :: mat(:,:) , b(:)
    REAL(KIND=DP), INTENT(INOUT) :: vec(:)
    REAL(KIND=DP), INTENT(OUT) :: bigdif     ! We'll use to determine error

    bigdif=0
    ! Apply step gauss-seidel
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

        ! Store in bigdif the biggest difference between components of vec between steps
        bigdif = MAX(bigdif, ABS(old - vec(i)))

    END DO
END SUBROUTINE pasitoGS
END MODULE gauss_seidel_pas

! Undefined number of steps. We'll just make it stop when it reaches the precision (tolerance) we want
MODULE gauss_seidel_tol
    USE gauss_seidel_pas
    IMPLICIT NONE (type, external)

    INTEGER :: numsteps
    REAL(KIND=DP) :: error

CONTAINS
SUBROUTINE gauseideltol(mat,vec,b,tolerance,stepslimit,print)
    REAL(KIND=DP), INTENT(IN) :: mat(:,:) , b(:) , tolerance
    INTEGER, INTENT(IN) :: stepslimit ! Just in case it doesn't converge, to make sure it stops at some point
    LOGICAL, INTENT(IN) :: print
    REAL(KIND=DP), INTENT(INOUT) :: vec(:) 

    ! Check whether sizes of matrix and vectors are compatible
    rows = SIZE(mat, 1)
    cols = SIZE(mat, 2)
    IF (SIZE(vec) /= rows) THEN
        PRINT *, "Error: Vector size must match the number of rows of the matrix."
        STOP
    ELSE IF (SIZE(b) /= rows) THEN
        PRINT *, "Error: Vector size must match the number of rows of the matrix."
        STOP
    END IF

    numsteps=0 ! We'll count how many steps are needed
    ! Apply Gauss-Seidel as long as error is bigger than tolerance
    DO WHILE (numsteps .LT. stepslimit)
        CALL pasitoGS(mat,vec,b,error)    ! Starting from given vector, we'll modify vector vec applying the method
        numsteps=numsteps+1

        IF (print) THEN
            IF (error .EQ. 0) THEN      ! this is the solution
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
