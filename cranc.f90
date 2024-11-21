Program cranc

   
    integer :: N=101
    integer :: i
    real :: Az 
    real :: At
    integer :: h,g,j
    integer :: c=490
    real :: Tc = 36.5
    real :: cv = 3686
    real :: r = 1081
    real :: k=0.56
    real :: s=0.472
    real :: v=40
    real :: L=0.02
    real :: p 
    real :: a
    real :: b
    real :: A_z=0.49
    real :: Tc_n
    REAL, ALLOCATABLE :: T(:,:)
    Az = 1.0/(N-1)
    At = A_z*(Az**2)
    p = (s*v**2)/(2*L**2)
    a=k/(cv*r)
    b=p/(r*cv)
    Tc_n=(Tc*a)/(b*L)
    c=int(0.025/At)
    allocate(T(0:N-1,0:N-1))

    do i = 1, N-1
         T(i,i)=-A_z/2
         T(i,i+1)=A_z/2
         T(i,i+2)=-A_z/2
    end do

    write(*, *) T(5,:)
end program cranc


