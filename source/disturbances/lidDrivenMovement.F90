    subroutine lidDrivenMovement(nx,ny,nz,x,y,z,t,u,U0)

    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real*8, dimension(nx), intent(in) :: x
    real*8, dimension(ny), intent(in) :: y
    real*8, dimension(nz), intent(in) :: z
    real*8, intent(in) :: t
    real*8, dimension(nx,ny,nz),intent(inout) :: u
    real*8, intent(in) :: U0
    
    integer :: i
    
    do i = 1,nx
        u(i,:,:) = U0 * (1- (2 * (x(i)-x(1)) / (x(nx)-x(1)) - 1)**18)**2
    enddo
    
    end subroutine
