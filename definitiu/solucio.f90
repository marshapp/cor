program solucio

implicit none(type, external) 
integer :: i, j, N_t
integer, parameter :: N_x = 101
real :: delta_t, alpha, lambda, P_ext, beta, delta_x
real, parameter :: L = 0.02, V = 40, c_v = 3686, rho = 1081, k = 0.56, sigma = 0.472, T_c=36.5, gamma=0.25, t_a=0.025
real, allocatable :: T(:,:) 

    delta_x = 1/(real(N_x)-1)

    alpha = k/(c_v*rho)

    P_ext = ((V**2)*sigma)/(2*(L**2))

    lambda = P_ext/(rho*c_v)    

    beta = T_c*((alpha)/(lambda*(L**2)))

    delta_t = gamma*(delta_x**2)

    N_t=int(t_a/delta_t)

    allocate(T(0:N_x-1, 0:N_t))

    !Apliquem la condició inicial:
    do j = 0, N_x-1
        T(j,0) = beta
    end  do

    !Apliquem les condicions de contorn:
    do i = 0, N_t
        T(0,i) = beta
        T(N_x-1,i) = beta
    end  do

    !Trobem els demés punts utilitzant la iteració del mètode:
    outer: do i = 0, N_t-1 !Excloem i=N_t per evitar sortir-nos de la matriu.

        do j = 1, N_x-2 !Excloem j=0 i j=N_x-1, per evitar sortir-nos de la matriu.
            T(j,i+1)=(T(j+1,i)-2*T(j,i)+T(j-1,i))*gamma+delta_t+T(j,i)
        end do

        do j = 0, N_x-1
            if (T(j,i+1)*((lambda*(L**2))/alpha) .ge. 80) then
                write(*,*) ((i+1)*delta_t)*((L**2)/alpha)
                exit outer
            end if
        end do

        do j = 0, 37
            if (T(j,i+1)*((lambda*(L**2))/alpha) .ge. 50) then
                write(*,*) ((i+1)*delta_t)*((L**2)/alpha)
                exit outer
            end if
        end do

        do j = 63, N_x-1
            if (T(j,i+1)*((lambda*(L**2))/alpha) .ge. 50) then
                write(*,*) ((i+1)*delta_t)*((L**2)/alpha)
                exit outer
            end if
        end do

    end do outer

end program solucio