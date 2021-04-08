    subroutine periodicSine(nx,ny,nz,x,y,z,t,u,omega,amplitude)

    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real*8, dimension(nx), intent(in) :: x
    real*8, dimension(ny), intent(in) :: y
    real*8, dimension(nz), intent(in) :: z
    real*8, intent(in) :: t
    real*8, dimension(nx,ny,nz),intent(inout) :: u
    real*8, intent(in) :: omega
    real*8, intent(in) :: amplitude
    
    integer :: i
    real :: alpha, pi
    
    pi = 4.d0*datan(1.d0)
    alpha = 2.d0*pi/(x(nx)-x(1))
    
    do i = 1,nx
        u(i,:,:) = u(i,:,:) + amplitude * dsin(alpha*(x(i)-x(1))) * dsin(omega*t)
    enddo
    
    end subroutine
