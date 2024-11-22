program jacobi
    integer, parameter :: N=3,len=1000
    real :: A(N,N), b(n), x(n), x_1(n)
    integer :: iter, i,j,k
    real :: sum, error
    real, parameter :: tolerancia=1.0e-9

        
    b= [2,3,-1] !escriu el vector del terme independent, b    
    A = reshape([10,-1,9, &
                 -1,5,4, &
                 1,2,7], [n, n])!definim la matriu A

    x= [0.0, 0.0, 0.0] !suposició inicial per x
    x_1= [0.0,0.0,0.0]
    
    print *, "ja estic resolent el teu sistema"
    
    do iter = 1, len
    
        x_1= x !la x de l'anterior iteració ara és x_1
        do i = 1, N !suma per cada x(i) 
            sum =0
            do j = 1, N
                if (j/=i) then
                    sum = sum + A(i,j)*x_1(j)
                end if
            end do 
            x(i)=(b(i)-sum)/A(i,i) !genera la nova x
        end do
        error = maxval(abs(x - x_1))

        if (error < tolerancia) then
            print*, "ha convergit en", iter, "iteracions"
            exit
        end if
        if (iter==len) then
            print*, "no ha convergit"
            stop 
        end if
        
    end do
        
    if (iter/=len) then
        write(*, *) "la solució al teu sistema és:"
        do i = 1, n
            print *, "x(", i, ") =", x(i)
        end do
    end if 

end program jacobi