program euler_explicit
    implicit none(type, external) 
    integer :: i, j, w, m, s_values(3) !u
    integer, parameter :: N = 101
    real :: T_0, delta_x, delta_t, gamma, alpha, beta, P_ext
    real, parameter :: L = 0.02, V = 40, c_v = 0.56, rho = 1081, k = 0.56, sigma = 0.472, T_c=36.5
    real :: c_values(3)
    real :: R(0:N-1,6)
    real, allocatable :: A(:,:)  ! Matriz dinámica
    !real :: A(0:N-1,0:490)  !estaria be dirll T a la matriu

    c_values=[0.51, 0.49, 0.25]

    delta_x = 1/(real(N)-1)
    
    do i = 1, 3
        delta_t = c_values(i)*(delta_x**2)
        s_values(i)=int(0.025/delta_t)+1
    end do

    delta_x = 1/(real(N)-1)

    alpha = k/(c_v*rho)

    P_ext = ((V**2)*sigma)/(2*(L**2))
    beta = P_ext/(rho*c_v)    


    T_0 = T_c*((alpha)/(beta*((L)**2)))

    OPEN(NEWUNIT=w, FILE='resultat.txt', STATUS='UNKNOWN', ACTION='WRITE')


    do m=1, 3

        allocate(A(0:N-1, 0:s_values(m)))
        A=0.0 !aleshores aixo es necessari???
        delta_t = c_values(m)*(delta_x**2)
        gamma = delta_t/(delta_x)**2

        do j = 0, N-1
            A(j,0) = T_0
        end  do

        do i = 0, s_values(m)
            A(0,i) = T_0
            A(N-1,i) = T_0
        end  do

        do i = 0, s_values(m)-1
            do j = 1, N-2
                A(j,i+1)=(A(j+1,i)-2*A(j,i)+A(j-1,i))*gamma+delta_t+A(j,i)
            end do
        end do

        do j = 0, N-1
            R(j,2*m-1) = (j*delta_x)*L
            R(j,2*m) = A(j, s_values(m))*((beta*(L**2))/alpha)
        end do

        deallocate(A)
 
    end do    

    do j = 0, N-1
        write(w, *) R(j,1), R(j,2), R(j,3), R(j,4), R(j,5), R(j,6)
    end do

    close(w)



end program euler_explicit



        

