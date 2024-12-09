Program crancgif

    implicit none(type, external)
    integer,parameter :: N=101 
    real, parameter :: Az =1.0/(N), Tc = 36.5, cv = 3686, r = 1081, k=0.56, s=0.472, v=40,  L=0.02, A_z=0.25
    real,parameter :: At= A_z*(Az**2)
    integer,parameter :: c=int(0.025/At), d=c/9
    real,parameter :: p = (s*v**2)/(2*L**2)
    real, parameter :: alpha = k/(cv*r)
    real,parameter :: beta=p/(r*cv)
    real,parameter :: Tc_n=(Tc*alpha)/(beta*L**2)
    
    
    integer, parameter :: len=1000
    real :: A(n-1,n-1), b(n-1), x(0:n,0:c), x_1(0:n,0:c)
    integer :: iter, i,j,t,h,g
    real :: sum, error
    real, parameter :: tolerancia=1.0e-9
    character(len=20) :: file_name
    !x és una matriu de temperatures, cada columna té les temperatures al llarg de l'espai per un temps
    

    
    do i=1,n-1 !inicialitzem la matriu A tota a 0
        do j=1,n-1
            A(i,j)= 0
        end do
    end do
    do i=1,n-1 ! inicialitzem la matriu A tota a 0
        do j=1,c
            x(i,j)= 0
        end do
    end do

    !construim la matriu A segons el mètode de crank-nicolson
    A(1,1)= A_z+1 
    A(1,2)= -A_z/2
    A(N-1,N-2)= -A_z/2
    A(N-1,N-1)=A_z+1
    
    do i = 2, N-2
        A(i,i-1)=-A_z/2
        A(i,i)=(A_z+1)
        A(i,i+1)=-A_z/2
    end do
    !Apliquem les condicions de contorn i les condicions inicials
    do h = 0,N !a temps 0 la temp és tc a tot l'espai
        x(h,0)= Tc_n 
        x_1(h,0)=Tc_n
    end do

    do h=0,c
        x(0,h)=Tc_n
        x(N,h)=Tc_n
        x_1(0,h)=Tc_n
        x_1(N,h)=Tc_n
    end do
    !resolem el sistema per el métode de jacobi
    do t= 0, c-1 !c-1 perque et donara x a tems t+1
        do j = 2,N-2
            !per cada pas temporal construim la matriu b, el terme independent
            B(1)=(A_z/2)*x(2,t)-(A_z-1)*x(1,t)+(A_z/2)*Tc_n+At+(A_z/2)*Tc_n
            B(N-1)=(A_z/2)*Tc_n-(A_z-1)*x(n-1,t)+(A_z/2)*x(n-2,t)+At+(A_z/2)*Tc_n
            B(j)=(A_z/2)*x(j+1,t)-(A_z-1)*x(j,t)+(A_z/2)*x(j-1,t)+At
        end do
    
        do iter = 1, len
        
            x_1(:,t+1)= x(:,t+1) !la x de l'anterior iteració ara és x_1
            do i = 1, N-1 !suma per cada x(i) 
                sum =0
                do j = 1, N-1
                    if (j/=i) then
                        sum = sum + A(i,j)*x_1(j,t+1)
                    end if
                end do 
                x(i,t+1)=(b(i)-sum)/A(i,i) !genera la nova x
            end do
            error = maxval(abs(x(:,t+1) - x_1(:,t+1)))

            if (error < tolerancia) then
                exit
            end if
            if (iter==len) then
                stop 
            end if
            
        end do
    end do   
        

        


    !Per cada temps generem un fitxer txt amb les dades de les temperatures a tot l'espai 
    do t = 0,d

        WRITE(file_name, '("d_", I0, ".txt")') t

        open(unit=t, file=file_name, status="replace")   
        do i = 0, n
            do j=0,n
                write(t, '(F10.6, 1X, F10.6, 1X, F10.6)') (i*Az)*l, (j*Az)*l, x(i,t*9)*((beta*(l**2))/alpha) 
            end do
            write(t,*)
        end do


        close(unit=t)

    end do

end program crancgif


 