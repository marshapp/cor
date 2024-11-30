
! kth step in Gauss-Seidel method ; **it's the same as in gauss_seidel_tolerance
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


! We say exactly how many iterations we'll want to run
MODULE gauss_seidel_ss
    USE gauss_seidel_pas
    IMPLICIT NONE (type, external)

    INTEGER :: numsteps
    REAL(KIND=DP) :: error

CONTAINS
SUBROUTINE gauseidelss(mat,vec,b,tolerance,maxsteps,print)
    REAL(KIND=DP), INTENT(IN) :: mat(:,:) , b(:) , tolerance
    INTEGER, INTENT(IN) :: maxsteps
    LOGICAL, INTENT(IN) :: print
    REAL(KIND=DP), INTENT(INOUT) :: vec(:)

    ! Check whether sizes of matrix and vectors are compatible
    rows = SIZE(mat, 1)
    cols = SIZE(mat, 2)
    IF (SIZE(vec) /= rows) THEN
        WRITE(*,*) "Error: Vector size must match the number of rows of the matrix."
        STOP
    ELSE IF (SIZE(b) /= rows) THEN
        WRITE(*,*) "Error: Vector size must match the number of rows of the matrix."
        STOP
    END IF

    numsteps=0 ! We'll keep track of how many iterations we've done
    ! Apply gauss-seidel as long as we haven't reached the maximum number of steps
    DO WHILE(numsteps .LT. maxsteps)
        CALL pasitoGS(mat,vec,b,error)    ! Starting from given vector, we'll modify vector vec applying the method
        numsteps=numsteps+1

        IF (error .EQ. 0) THEN   ! We found the solution
            EXIT
        END IF
    END DO

    IF (print) THEN
        IF (error .EQ. 0) THEN
            WRITE(*,'(A,1X,I0,1X,A)') 'The method converged in ', numsteps-1 , ' iterations. The solution is: '
        ELSE IF (error .LT. tolerance) THEN     ! We've reached the precision we want 
            WRITE(*,*) 'The method converges; the desired precision was reached. The approximate solution is: '
        END IF

        DO i=1,rows
            WRITE(*,*) '$x_{',i,'} = $', vec(i)
        END DO
    END IF

    IF (error .GT. tolerance) THEN
        WRITE(*,*) 'The method does not converge.'
    END IF

END SUBROUTINE gauseidelss
END MODULE gauss_seidel_ss