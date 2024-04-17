    subroutine turbulenceGenerator(nx,ny,nz,x,y,z,u,v,w,du,dv,dw,epsilon)

    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real*8, dimension(nx), intent(in) :: x
    real*8, dimension(ny), intent(in) :: y
    real*8, dimension(nz), intent(in) :: z
    real*8, dimension(nx,ny,nz),intent(in) :: u, v, w
    real*8, dimension(nx,ny,nz),intent(inout) :: du, dv, dw
    real*8, intent(in) :: epsilon
    
    du = du + epsilon*u
    dv = dv + epsilon*v
    dw = dw + epsilon*w
    
    end subroutine
