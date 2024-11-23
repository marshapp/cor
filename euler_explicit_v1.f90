program euler_explicit

    implicit none(type, external) 
    integer :: i, j, m, u, h(3)
    integer, parameter :: N_x = 101
    real :: T_0, delta_x, delta_t, alpha, lambda, P_ext, beta, gamma(3), R(0:N_x-1,6)
    real, parameter :: L = 0.02, V = 40, c_v = 0.56, rho = 1081, k = 0.56, sigma = 0.472, T_c=36.5
    real, allocatable :: A(:,:) 

    gamma=[0.51, 0.49, 0.25]

    delta_x = 1/(real(N_x)-1)
    
    do i = 1, 3
        delta_t = gamma(i)*(delta_x**2)
        h(i)=int(0.025/delta_t)+1
    end do

    alpha = k/(c_v*rho)

    P_ext = ((V**2)*sigma)/(2*(L**2))

    lambda = P_ext/(rho*c_v)    

    beta = T_c*((alpha)/(lambda*(L**2)))

    T_0 = T_c*((alpha)/(lambda*((L)**2)))

    open(newunit=u, file='resultats_euler_explicit.txt', status='unknown', action='write')

    do m=1, 3

        allocate(A(0:N_x-1, 0:h(m)))
    
        delta_t = gamma(m)*(delta_x**2)

        do j = 0, N_x-1
            A(j,0) = T_0
        end  do

        do i = 0, h(m)
            A(0,i) = T_0
            A(N_x-1,i) = T_0
        end  do

        do i = 0, h(m)-1
            do j = 1, N_x-2
                A(j,i+1)=(A(j+1,i)-2*A(j,i)+A(j-1,i))*gamma(m)+delta_t+A(j,i)
            end do
        end do

        do j = 0, N_x-1
            R(j,2*m-1) = (j*delta_x)*L
            R(j,2*m) = A(j, h(m))*((lambda*(L**2))/alpha)
        end do

        deallocate(A)
 
    end do    

    do j = 0, N_x-1
        write(u, *) R(j,1), R(j,2), R(j,3), R(j,4), R(j,5), R(j,6)
    end do

    close(u)

end program euler_explicit