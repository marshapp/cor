program explicit

    implicit none(type, external) 
    integer :: i, j, m, u, w, N_t
    integer, parameter :: N_z = 101
    real :: dz, dt, alpha, lambda, P_ext, beta, gamma(3), R(0:N_z-1,3)
    real, parameter :: L = 0.02, V = 40, c_v = 3686, rho = 1081, k = 0.56, sigma = 0.472, T_c=36.5, t_a=0.025
    real, allocatable :: T(:,:) 
    
    dz = 1/(real(N_z)-1) !Definim l'interval espacial

    !Definim els diferents valors de gamma o mallats que utilitzarem:
    gamma=[0.51, 0.49, 0.25]
    
    !Definim diferents constants del problema:
    alpha = k/(c_v*rho)

    P_ext = ((V**2)*sigma)/(2*(L**2))

    lambda = P_ext/(rho*c_v)    

    beta = T_c*((alpha)/(lambda*(L**2)))


    do m=1, 3 !Fem un bucle per tenir en compte els 3 valors de gamma
        
        dt = gamma(m)*(dz**2) !Definim el dt associat a aquest valor de gamma

        N_t=int(t_a/dt) !Definim el nombre de passos temporals que depèn de dt i, per tant, també de gamma

        allocate(T(0:N_z-1, 0:N_t)) !Definim les dimensions de la matriu que emmagatzema 
        !els valors de la temperatura a cada temps i posició per aquest valor de gamma

        !Apliquem la condició inicial:
        do j = 0, N_z-1
            T(j,0) = beta
        end  do

        !Apliquem les condicions de contorn:
        do i = 0, N_t
            T(0,i) = beta
            T(N_z-1,i) = beta
        end  do

        !Trobem els demés punts utilitzant l'equació corresponent pel temps i+1 i la posició j:
        do i = 0, N_t-1 !Excloem i=N_t per evitar sortir-nos de la matriu
            do j = 1, N_z-2 !Excloem j=0 i j=N_x-1 per evitar sortir-nos de la matriu
                T(j,i+1)=(T(j+1,i)-2*T(j,i)+T(j-1,i))*gamma(m)+dt+T(j,i)
            end do
        end do

        !Guardem els resultats sense normalitzar en una nova matriu per després facilitar el tractament dels mateixos:
        do j = 0, N_z-1
            R(j,m) = T(j, N_t)*((lambda*(L**2))/alpha)
        end do

        deallocate(T) !Reiniciem la matriu per a emmagatzemar els valors de la temperatura pel següent valor de gamma
 
    end do    


    open(newunit=u, file='explicit.dat', status='unknown', action='write')

    !Escribim en un fitxer els punts espacials sense normalitzar i els resultats corresponents a cada gamma, 
    !en 4 columnes diferents per facilitar el seu graficat:
    do j = 0, N_z-1
        write(u, *) (j*dz)*L, R(j,1), R(j,2), R(j,3), f(200, j*dz)
    end do

    close(u)


    open(NEWUNIT=w, file='error_explicit.dat', status='unknown', action='write')

    !Escribim en un fitxer els punts espacials sense normalitzar i els errors absoluts corresponents a gamma=0,49
    !i a gamma=0,25 en 3 columnes diferents per facilitar el seu graficat:
    do j = 0, N_z-1
        write(w, *) (j*dz)*L, abs(R(j,2)-f(200, j*dz)), abs(R(j,3)-f(200, j*dz))
    end do

    close(w)

    contains 
    !Aquí definim la funció corresponent a la solució analítica pel temps normalitzat t_a
    !que pren per arguments la posició normalitzada i el nombre de termes del sumatori que apareix a la seva expressió    
        real function f(N, z_norm)

            implicit none(type, external)
            integer, intent(in) :: N
            real, intent(in) :: z_norm
            real, parameter :: pi=2*acos(0.0)
            integer :: k
            real :: sumatori, f_norm

            sumatori=0 

            do k = 1, N
                sumatori = sumatori + ((1-exp(-((2*k-1)**2)*(pi**2)*t_a))/((2*k-1)**3))*sin((2*k-1)*pi*z_norm)
            end do

            f_norm=beta+(4/(pi**3))*sumatori
            f=(f_norm*lambda*(L**2))/alpha !Desnormalitzem el resultat de sortida

        end function

end program explicit
