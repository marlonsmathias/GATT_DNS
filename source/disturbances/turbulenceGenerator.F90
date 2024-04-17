    subroutine turbulenceGenerator(nx,ny,nz,x,y,z,u,v,w,du,dv,dw,epsilon)

    use mpi

    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real*8, dimension(nx), intent(in) :: x
    real*8, dimension(ny), intent(in) :: y
    real*8, dimension(nz), intent(in) :: z
    real*8, dimension(nx,ny,nz),intent(in) :: u, v, w
    real*8, dimension(nx,ny,nz),intent(inout) :: du, dv, dw
    real*8, intent(in) :: epsilon

    real*8 :: E_local, E
    integer :: N_local, N
    integer :: ierror
    
    N_local = nx * ny * nz
    E_local = sum(U*U + V*V + W*W)

    call MPI_ALLREDUCE(E_local, E, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierror)
    call MPI_ALLREDUCE(N_local, N, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierror)

    E = E/N

    du = du + 0.5*epsilon/E*u
    dv = dv + 0.5*epsilon/E*v
    dw = dw + 0.5*epsilon/E*w
    
    end subroutine
