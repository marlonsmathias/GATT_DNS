    ! This routine is used as boundary condition to hold the flow at its initial state for U,V,W,R and E, for example when a blasius profile inflow must be used
	
	subroutine holdInlet(nx,ny,nz,x,y,z,t,u,v,r,w,e,holdDensityFloat)

    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real*8, dimension(nx), intent(in) :: x
    real*8, dimension(ny), intent(in) :: y
    real*8, dimension(nz), intent(in) :: z
    real*8, intent(in) :: t, holdDensityFloat
    real*8, dimension(nx,ny,nz),intent(inout) :: u,v,r,w,e
    
    logical, save :: firstCall = .TRUE.
	logical, save :: holdDensity
    real*8, dimension(:,:,:), allocatable, save :: uSave,vSave,rSave,wSave,eSave
	
    if (firstCall) then
		firstCall = .FALSE.
		
		allocate(uSave(nx,ny,nz))
		allocate(vSave(nx,ny,nz))
		allocate(wSave(nx,ny,nz))
		allocate(eSave(nx,ny,nz))
		
		uSave = u
		vSave = v
		wSave = w
		eSave = e
		
		if (holdDensityFloat.gt.0.5d0) then
			holdDensity = .TRUE.
			allocate(rSave(nx,ny,nz))
			rSave = r
		else
			holdDensity = .FALSE.
		endif
		
	else
		u = uSave
		v = vSave
		w = wSave
		e = eSave
		
		if (holdDensity) then
			r = rSave
		endif
	endif
    
    end subroutine
