    subroutine turbulenceGenerator(nx,ny,nz,x,y,z,u,v,w,e,du,dv,dw,de,epsilon,e0)

    use mpi

    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real*8, dimension(nx), intent(in) :: x
    real*8, dimension(ny), intent(in) :: y
    real*8, dimension(nz), intent(in) :: z
    real*8, dimension(nx,ny,nz),intent(in) :: u, v, w, e
    real*8, dimension(nx,ny,nz),intent(inout) :: du, dv, dw, de
    real*8, intent(in) :: epsilon

    real*8 :: ET_local, ET_sum, U_local, U_sum, V_local, V_sum, W_local, W_sum, E_local, E_sum
    integer :: N_local, N
    integer :: ierror
    real*8 :: gamma, e0
    
    N_local = nx * ny * nz
    ET_local = sum(U*U + V*V + W*W)
    U_local = sum(U)
    V_local = sum(V)
    W_local = sum(W)
    E_local = sum(E)

    gamma = 10*epsilon

    call MPI_ALLREDUCE(ET_local, ET_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierror)

    call MPI_ALLREDUCE(U_local, U_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierror)
    call MPI_ALLREDUCE(V_local, V_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierror)
    call MPI_ALLREDUCE(W_local, W_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierror)
    call MPI_ALLREDUCE(E_local, E_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierror)

    call MPI_ALLREDUCE(N_local, N, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierror)

    du = du + 0.25*epsilon*N/ET_sum*u - gamma*U_sum/N
    dv = dv + 0.25*epsilon*N/ET_sum*v - gamma*V_sum/N
    dw = dw + 0.25*epsilon*N/ET_sum*w - gamma*W_sum/N

    de = de - gamma*(E_sum/N-e0)
    
    end subroutine
