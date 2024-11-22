program jacobi
    integer, parameter :: N=3,len=20
    real :: A(N,N), b(n), x(n), x_1(n)
    integer :: iter, i,j,k
    real :: sum, error
    real, parameter :: tolerancia=1.0e-9

        
    b= [15,10,10] !escriu el vector del terme independent, b    
    A = reshape([4.0, -1.0, 0.0, &
             -1.0, 4.0, -1.0, &
              0.0, -1.0, 4.0], [n, n])!definim la matriu A

    x= reshape([0,0,0], [n]) !suposició inicial per x
    x_1= reshape([0,0,0], [n])
    
    print *, "ja estic resulent el teu sistema"
    
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
    end do
        
        
    write(*, *) "la solució al teu sistema és:"
    do i = 1, n
        print *, "x(", i, ") =", x(i)
    end do


end program jacobi