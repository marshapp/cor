program errcranc2

        implicit none(type, external)

        integer :: w
        integer, parameter :: n=101
        real, parameter :: Az2 =1.0/(N), l=0.02, leng = 1000
        real, parameter :: gamma(2)

        gamma = [1.0, 0.5]

    
        open(unit=13, file="cranc_error.txt", status="replace")   
        do w = 0, N
            !escriu: Az, Az/l, Temp(Az,t=c) funcio analitica(N=1000), Temp(Az,t=c) cranc, error absolut
            write(13, *) (w*Az2)*l, (w*Az2), abs(f(leng, w*Az2)-cranc(a))
            !, abs(f(len, w*Az)-x(w,c)*((beta*(l**2))/alpha))
            !f(len, w*Az), x(w,c)*((beta*(l**2))/alpha)
            
        end do
    
        close(13)
    
    contains

    real function f(m, x_norm)
        implicit none(type, external)
        integer, intent(in) :: m
        real, intent(in) :: x_norm
        real, parameter :: t_norm=0.025, pi=2*acos(0.0)
        real :: sumatori, f_norm
        real, parameter :: n=101
        real, parameter :: Az =1.0/(N), Tc = 36.5, cv = 3686, r = 1081, k=0.56, s=0.472, v=40,  L=0.02 
        real,parameter :: p = (s*v**2)/(2*L**2)
        real, parameter :: alpha = k/(cv*r)
        real,parameter :: beta = p/(r*cv)
        real, parameter :: betta = Tc*((alpha)/(beta*(L**2)))
        integer :: o
    
        sumatori=0 
        do o = 1, m
            sumatori = sumatori + ((1-exp(-((2*k-1)**2)*(pi**2)*t_norm))/((2*k-1)**3))*sin((2*k-1)*pi*x_norm)
        end do
    
        f_norm=betta+(4/(pi**3))*sumatori
        f=(f_norm*beta*(L**2))/alpha
    
    end function
    
    real function cranc(A_z,x_normc)

        implicit none(type, external)
        real, intent(in) :: A_z
        real, intent(in) :: x_normc
        integer,parameter :: N=101 !al començar a n=0 hi haura n+1 punts
        integer, parameter :: N_x = 101
        real, parameter :: Az =1.0/(N), Tc = 36.5, cv = 3686, r = 1081, k=0.56, s=0.472, v=40,  L=0.02 
        real :: At= A_z*(Az**2)
        integer :: c = int(0.025/At)
        real,parameter :: p = (s*v**2)/(2*L**2)
        real, parameter :: alpha = k/(cv*r) 
        real,parameter :: beta = p/(r*cv)
        !real, parameter :: betta = Tc*((alpha)/(beta*(L**2)))
        real,parameter :: Tc_n=(Tc*alpha)/(beta*L**2)
        real :: cranc_norm
        
        integer, parameter :: len=1000
        real :: A(n-1,n-1), b(n-1), x(0:n,0:c), x_1(0:n,0:c)
        integer :: iter, i, j, t, h, w
        real :: sum, error
        real, parameter :: delta_x = 1/(real(N_x)-1)
        real, parameter :: tolerancia=1.0e-9
        
        do i=1,n-1 !matriu a  tot a 0, matriu x tot a 0 (després: menys a t=0  que val Tc)
            do j=1,n-1
                A(i,j)= 0
            end do
        end do
        do i=1,n-1 ! x tot a 0 (després: menys a t=0  que val Tc)
            do j=1,c
                x(i,j)= 0
            end do
        end do
    
        A(1,1)= A_z+1 
        A(1,2)= -A_z/2
        A(N-1,N-2)= -A_z/2
        A(N-1,N-1)=A_z+1
    
        do i = 2, N-2
            A(i,i-1)=-A_z/2
            A(i,i)=(A_z+1)
            A(i,i+1)=-A_z/2
        end do
    
        do h = 0,N !a temps 0 la temp és tc a tot l'espai
            x(h,0)= Tc_n !ci
            x_1(h,0)=Tc_n
        end do
    
        do h=0,c
            x(0,h)=Tc_n
            x(N,h)=Tc_n
            x_1(0,h)=Tc_n
            x_1(N,h)=Tc_n
        end do
    
        do t= 0, c-1 !c-1 perque et donara x a tems t+1
            do j = 2,N-2
                B(1)=(A_z/2)*x(2,t)-(A_z-1)*x(1,t)+(A_z/2)*Tc_n+At+(A_z/2)*Tc_n
                B(N-1)=(A_z/2)*Tc_n-(A_z-1)*x(n-1,t)+(A_z/2)*x(n-2,t)+At+(A_z/2)*Tc_n
                B(j)=(A_z/2)*x(j+1,t)-(A_z-1)*x(j,t)+(A_z/2)*x(j-1,t)+At
            end do
            !print *, "ja estic resolent el teu sistema"
        
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
                    !print*, "ha convergit en", iter, "iteracions"
                    exit
                end if
                if (iter==len) then
                    print*, "no ha convergit ): sap greu tant jove"
                    stop 
                end if
                
            end do
        end do   

        cranc_norm = x(i,c)
        cranc = cranc_norm*((beta*(l**2))/alpha)
        
    end function

    end program errcranc2