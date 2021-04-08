    subroutine packet_2d(nx,ny,nz,x,y,z,t,u,omega,nmodes,amplitude)

    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real*8, dimension(nx), intent(in) :: x
    real*8, dimension(ny), intent(in) :: y
    real*8, dimension(nz), intent(in) :: z
    real*8, intent(in) :: t
    real*8, dimension(nx,ny,nz),intent(inout) :: u
    real*8, intent(in) :: omega    ! fundamental radial frequency / spanwise frequency
    real*8, intent(in) :: amplitude ! of each mode
    real*8, intent(in) :: nmodes  ! number of Fourier modes as a real number for compatibilty with runDisturbance
    
    integer :: ix,iz,imode,nmode,jmode,bmode
    real :: pi,x0,Tp,Tend,alpha,sigmaT,sigmaX
    complex :: iImag,aux
    
    iImag = (0.0,1.0)
    pi = 4.d0*datan(1.d0)
    alpha = 2.d0*pi/32d0;
    nmode=int(nmodes)
    x0=(x(nx)+x(1))/2;
    Tend = (2.d0*pi/omega)*2.5d-1
    Tp= Tend*5d-1;
    sigmaT = Tend*2.5d-1;
    sigmaX = (x(nx)-x(1))/8.d0;
    if(t.gt.Tend) then
      return
    else
      do ix=1,nx
        do iz=1,nz
          do imode=1,nmode
			  aux = exp(iImag*(imode*omega*(t-Tp)))
              u(ix,:,iz) = u(ix,:,iz) + dexp(-(t-Tp)**(2.d0)/(2.d0*sigmaT**(2.d0)))*dexp(-(x(ix)-x0)**(2.d0)/(2.d0*sigmaX**(2.d0)))*real(aux)/(nmodes)
          enddo !imode=1,nmode
          u(ix,:,iz)= amplitude*u(ix,:,iz);
          enddo !iz=1,nz
      enddo !ix=1,nx
    end if
    end subroutine
