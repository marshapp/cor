program solucio

    implicit none(type, external) 
    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
    
    integer :: i, j, N_t
    integer, parameter :: N_z = 101
    real(KIND=DP) :: dz, dt, alpha, lambda, P_ext, beta
    real(KIND=DP), parameter :: L = 0.02, V = 40, c_v = 3686, rho = 1081, k = 0.56, sigma = 0.472, T_c=36.5, gamma=0.25, t_a=0.025
    real(KIND=DP), allocatable :: T(:,:) 

    dz = 1/(real(N_z)-1) !Definim l'interval espacial

    !Definim diferents constants del problema:
    alpha = k/(c_v*rho)

    P_ext = ((V**2)*sigma)/(2*(L**2))

    lambda = P_ext/(rho*c_v)    

    beta = T_c*((alpha)/(lambda*(L**2)))

    dt = gamma*(dz**2) !Definim el dt associat a aquest gamma

    N_t=int(t_a/dt) !Definim el nombre de passos temporals que depèn de dt i, per tant, també de gamma

    allocate(T(0:N_z-1, 0:N_t)) !Definim les dimensions de la matriu que emmagatzema 
    !els valors de la temperatura a cada temps i posició

    !Apliquem la condició inicial:
    do j = 0, N_z-1
        T(j,0) = beta
    end  do

    !Apliquem les condicions de contorn:
    do i = 0, N_t
        T(0,i) = beta
        T(N_z-1,i) = beta
    end  do

    !Trobem els demés punts utilitzant la iteració del mètode:
    outer: do i = 0, N_t-1 !Excloem i=N_t per evitar sortir-nos de la matriu.

        do j = 1, N_z-2 !Excloem j=0 i j=N_z-1, per evitar sortir-nos de la matriu.
            T(j,i+1)=(T(j+1,i)-2*T(j,i)+T(j-1,i))*gamma+dt+T(j,i)
        end do

        !Comprovem si es compleix la condició de que la temperatura a ninguna part del teixit pot ser major o igual a 80 graus:
        do j = 0, N_z-1
            if (T(j,i+1)*((lambda*(L**2))/alpha) .ge. 80) then
                write(*,*) ((i+1)*dt)*((L**2)/alpha)
                exit outer
            end if
        end do

        !Comprovem si es compleix la condició de que la temperatura a la part sana del teixit no pot ser major o igual a 50 graus
        !en la part del teixit que va dels 0 cm als 0,75 cm
        !Si la temperatura supera o és igual a 80 graus en aquesta regió, s'imprimeix el temps desnormalitzat a la terminal i finalitza el programa
        do j = 0, 37
            if (T(j,i+1)*((lambda*(L**2))/alpha) .ge. 50) then
                write(*,*) ((i+1)*dt)*((L**2)/alpha)
                exit outer
            end if
        end do

        !Comprovem si es compleix la condició de que la temperatura a la part sana del teixit no pot ser major o igual a 50 graus
        !en la part del teixit que va dels 1,25 cm als 2 cm
        !Si la temperatura supera o és igual a 50 graus en aquesta regió, s'imprimeix el temps desnormalitzat a la terminal i finalitza el programa
        do j = 63, N_z-1
            if (T(j,i+1)*((lambda*(L**2))/alpha) .ge. 50) then
                write(*,*) ((i+1)*dt)*((L**2)/alpha) 
                exit outer
            end if
        end do

    end do outer

end program solucio
