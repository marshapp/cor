program prova

    implicit none(type, external)
    integer :: i, j, u
    integer, parameter :: N_x = 101
    real :: beta, alpha, lambda, P_ext, delta_x, xs(101)
    real, parameter :: L = 0.02, V = 40, c_v = 0.56, rho = 1081, k = 0.56, sigma = 0.472, T_c=36.5
    
    delta_x = 1/(real(N_x)-1)
    
    alpha = k/(c_v*rho)
    
    P_ext = ((V**2)*sigma)/(2*(L**2))
    
    lambda = P_ext/(rho*c_v) 
    
    beta = T_c*((alpha)/(lambda*(L**2)))
    
    open(newunit=u, file='resultats_analitica.txt', status='UNKNOWN', action='WRITE')
    
    do i = 0, 100
        xs(i+1)=i*delta_x
    end do
    
    do j = 1, N_x
        write(u, *) xs(j)*L, f(10, xs(j)), f(100, xs(j)), f(1000, xs(j))
    end do
    
    close(u)
    
    contains
    
    real function f(N, x_norm)
        implicit none(type, external)
        integer, intent(in) :: N
        real, intent(in) :: x_norm
        real, parameter :: t_norm=0.025, pi=2*acos(0.0)
        real :: sumatori, f_norm
        integer :: k
    
        sumatori=0 
        do k = 1, N
            sumatori = sumatori + ((1-exp(-((2*k-1)**2)*(pi**2)*t_norm))/((2*k-1)**3))*sin((2*k-1)*pi*x_norm)
        end do
    
        f_norm=beta+(4/(pi**3))*sumatori
        f=(f_norm*lambda*(L**2))/alpha
    
    end function
    
end program prova