program gif

    implicit none(type, external) 
    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
    
    integer :: i, j, N_t, temps, d
    integer, parameter :: N_z = 101
    real(KIND=DP) :: dz, dt, alpha, lambda, P_ext, beta
    real(KIND=DP), parameter :: L = 0.02, V = 40, c_v = 3686, rho = 1081, k = 0.56, sigma = 0.472, T_c=36.5, gamma=0.25, t_a=0.025
    real(KIND=DP), allocatable :: T(:,:) 
    character(len=20) :: file_name

    dz = 1/(real(N_z)-1) !Definim l'interval espacial

    !Definim diferents constants del problema:
    alpha = k/(c_v*rho)

    P_ext = ((V**2)*sigma)/(2*(L**2))

    lambda = P_ext/(rho*c_v)    

    beta = T_c*((alpha)/(lambda*(L**2)))

    dt = gamma*(dz**2) !Definim el dt associat a aquest gamma
    
    N_t=int(t_a/dt) !Definim el nombre de passos temporals que depèn de dt i, per tant, també de gamma

    d=N_t/9

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
    do i = 0, N_t-1 !Excloem i=N_t per evitar sortir-nos de la matriu.

        do j = 1, N_z-2 !Excloem j=0 i j=N_z-1, per evitar sortir-nos de la matriu.
            T(j,i+1)=(T(j+1,i)-2*T(j,i)+T(j-1,i))*gamma+dt+T(j,i)
        end do

    end do

    !Per cada temps generem un fitxer txt amb les dades de les temperatures a tot l'espai 
    do temps = 0, d

        WRITE(file_name, '("d_", I0, ".txt")') temps

        open(unit=temps, file=file_name, status="replace") 

        do i = 0, N_z-1
            do j=0, N_z-1
                write(temps, '(F10.6, 1X, F10.6, 1X, F10.6)') (i*dz)*L, (j*dz)*L, T(i,temps*9)*((lambda*(L**2))/alpha) 
            end do
            write(temps,*)
        end do

        close(unit=temps)

    end do

end program gif
