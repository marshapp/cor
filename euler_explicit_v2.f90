program euler_explicit

    implicit none(type, external) 
    integer :: i, j, m, u, w, N_t
    integer, parameter :: N_x = 101
    real :: delta_t, alpha, lambda, P_ext, beta, gamma(3), R(0:N_x-1,3), delta_x
    real, parameter :: L = 0.02, V = 40, c_v = 3686, rho = 1081, k = 0.56, sigma = 0.472, T_c=36.5, t_a=0.025
    real, allocatable :: T(:,:) 

    delta_x = 1/(real(N_x)-1)

    gamma=[0.51, 0.49, 0.25]

    alpha = k/(c_v*rho)

    P_ext = ((V**2)*sigma)/(2*(L**2))

    lambda = P_ext/(rho*c_v)    

    beta = T_c*((alpha)/(lambda*(L**2)))


    do m=1, 3

        delta_t = gamma(m)*(delta_x**2)

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
        do i = 0, N_t-1 !Excloem i=N_t per evitar sortir-nos de la matriu.
            do j = 1, N_x-2 !Excloem j=0 i j=N_x-1 per evitar sortir-nos de la matriu.
                T(j,i+1)=(T(j+1,i)-2*T(j,i)+T(j-1,i))*gamma(m)+delta_t+T(j,i)
            end do
        end do

        !Guardem els resultats sense normalitzar en una nova matriu per després facilitar el tractament dels mateixos:
        do j = 0, N_x-1
            R(j,m) = T(j, N_t)*((lambda*(L**2))/alpha)
        end do



        deallocate(T)
 
    end do    

    open(newunit=u, file='resultats_euler_explicit_v2.txt', status='unknown', action='write')

    !Escribim en un fitxer els punts espacials sense normalitzar i els resultats corresponents a cada gamma, 
    !en 4 columnes diferents per facilitar el seu graficat:
    do j = 0, N_x-1
        write(u, *) (j*delta_x)*L, R(j,1), R(j,2), R(j,3)
    end do

    close(u)


    open(NEWUNIT=w, file='resultats_error_v2.txt', status='unknown', action='write')

    !Escribim en un fitxer els punts espacials sense normalitzar i els errors absoluts corresponents a gamma=0.49
    !i a gamma=0.25 en 3 columnes diferents per facilitar el seu graficat:
    do j = 0, N_x-1
        write(w, *) (j*delta_x)*L, abs(R(j,2)-f(200, j*delta_x)), abs(R(j,3)-f(200, j*delta_x))
    end do

    close(w)

    contains

        real function f(N, x_norm)

            implicit none(type, external)
            integer, intent(in) :: N
            real, intent(in) :: x_norm
            real, parameter :: t_norm=0.025, pi=2*acos(0.0)
            integer :: k
            real :: sumatori, f_norm

            sumatori=0 

            do k = 1, N
                sumatori = sumatori + ((1-exp(-((2*k-1)**2)*(pi**2)*t_norm))/((2*k-1)**3))*sin((2*k-1)*pi*x_norm)
            end do

            f_norm=beta+(4/(pi**3))*sumatori
            f=(f_norm*lambda*(L**2))/alpha

        end function

end program euler_explicit