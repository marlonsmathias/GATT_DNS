    subroutine pressureGradient(nx,ny,nz,x,y,z,r,u,dr,du,dpdx)

    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real*8, intent(in) :: dpdx
    real*8, dimension(nx), intent(in) :: x
    real*8, dimension(ny), intent(in) :: y
    real*8, dimension(nz), intent(in) :: z
    real*8, dimension(nx,ny,nz),intent(in) :: r, u
    real*8, dimension(nx,ny,nz),intent(inout) :: dr, du

    du = du + dpdx/r
    
    end subroutine
