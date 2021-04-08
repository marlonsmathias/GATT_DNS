    module derivs
    
    use decomp_2d
    
    contains
    
    ! COMPUTE THE DERIVATIVES
    
    subroutine calcDerivsX(F,Fx,nDerivBlocksX,derivBlocksX, &
                        & derivnRHSx, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsDX, &
                        & periodicX)
                        
        ! INPUTS
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: F
        integer, intent(in)                   :: nDerivBlocksX
        integer, dimension(:,:), intent(in)   :: derivBlocksX
        integer, intent(in)                   :: derivnRHSx
        real*8,  dimension(:,:), intent(in)   :: derivsAX, derivsBX, derivsCX
        real*8,  dimension(:,:), intent(in)   :: derivsDX
        real*8,  dimension(:,:,:), intent(in) :: derivsRX
        integer, intent(in)                   :: periodicX
        
        ! OUTPUT
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Fx
        
        select case(periodicX)
        case(0)
            call solveSystemX(F,Fx,nDerivBlocksX,derivBlocksX,derivsAX,derivsBX,derivsCX,derivnRHSx,derivsRX)
        case(1)
            call solveSystemXper(F,Fx,nDerivBlocksX,derivBlocksX,derivsAX,derivsBX,derivsCX,derivsDX,derivnRHSx,derivsRX)
        end select
                        
    end subroutine
    
    subroutine calcDerivsY(F,Fy,nDerivBlocksY,derivBlocksY, &
                        & derivnRHSy, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsDY, &
                        & periodicY)
                        
        ! INPUTS
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: F
        integer, intent(in)                   :: nDerivBlocksY
        integer, dimension(:,:), intent(in)   :: derivBlocksY
        integer, intent(in)                   :: derivnRHSy
        real*8,  dimension(:,:), intent(in)   :: derivsAY, derivsBY, derivsCY
        real*8,  dimension(:,:), intent(in)   :: derivsDY
        real*8,  dimension(:,:,:), intent(in) :: derivsRY
        integer, intent(in)                   :: periodicY
        
        ! OUTPUT
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Fy
        
        ! INTERNALS
        real*8,  dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: F2, Fy2
        
        call transpose_x_to_y(F,F2)
        
        select case(periodicY)
        case(0)
            call solveSystemY(F2,Fy2,nDerivBlocksY,derivBlocksY,derivsAY,derivsBY,derivsCY,derivnRHSy,derivsRY)
        case(1)
            call solveSystemYper(F2,Fy2,nDerivBlocksY,derivBlocksY,derivsAY,derivsBY,derivsCY,derivsDY,derivnRHSy,derivsRY)
        end select
        
        call transpose_y_to_x(Fy2,Fy)
                        
                        
    end subroutine
    
    subroutine calcDerivsZ(F,Fz,nDerivBlocksZ,derivBlocksZ, &
                        & derivnRHSz, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDZ, &
                        & periodicZ)
                        
        ! INPUTS
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: F
        integer, intent(in)                   :: nDerivBlocksZ
        integer, dimension(:,:), intent(in)   :: derivBlocksZ
        integer, intent(in)                   :: derivnRHSz
        real*8,  dimension(:,:), intent(in)   :: derivsAZ, derivsBZ, derivsCZ
        real*8,  dimension(:,:), intent(in)   :: derivsDZ
        real*8,  dimension(:,:,:), intent(in) :: derivsRZ
        integer, intent(in)                   :: periodicZ
        
        ! OUTPUT
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Fz
        
        ! INTERNALS
        real*8,  dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: F2, Fz2
        real*8,  dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: F3, Fz3
        
        call transpose_x_to_y(F,F2)
        call transpose_y_to_z(F2,F3)
        
        select case(periodicZ)
        case(0)
            call solveSystemZ(F3,Fz3,nDerivBlocksZ,derivBlocksZ,derivsAZ,derivsBZ,derivsCZ,derivnRHSz,derivsRZ)
        case(1)
            call solveSystemZper(F3,Fz3,nDerivBlocksZ,derivBlocksZ,derivsAZ,derivsBZ,derivsCZ,derivsDZ,derivnRHSz,derivsRZ)
        end select
        
        call transpose_z_to_y(Fz3,Fz2)
        call transpose_y_to_x(Fz2,Fz)
                        
                        
    end subroutine
    
    subroutine calcDerivsXY(F,Fx,Fy,nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                        & derivnRHSx,derivnRHSy, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsDX, derivsDY, &
                        & periodicX, periodicY)
                        
        ! INPUTS
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: F
        integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksY
        integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksY
        integer, intent(in)                   :: derivnRHSx, derivnRHSy
        real*8,  dimension(:,:), intent(in)   :: derivsAX, derivsBX, derivsCX, derivsAY, derivsBY, derivsCY
        real*8,  dimension(:,:), intent(in)   :: derivsDX, derivsDY
        real*8,  dimension(:,:,:), intent(in) :: derivsRX, derivsRY
        integer, intent(in)                   :: periodicX, periodicY
        
        ! OUTPUT
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Fx, Fy
        
        ! INTERNALS
        real*8,  dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: F2, Fy2
        
        select case(periodicX)
        case(0)
            call solveSystemX(F,Fx,nDerivBlocksX,derivBlocksX,derivsAX,derivsBX,derivsCX,derivnRHSx,derivsRX)
        case(1)
            call solveSystemXper(F,Fx,nDerivBlocksX,derivBlocksX,derivsAX,derivsBX,derivsCX,derivsDX,derivnRHSx,derivsRX)
        end select
        
        call transpose_x_to_y(F,F2)
        select case(periodicY)
        case(0)
            call solveSystemY(F2,Fy2,nDerivBlocksY,derivBlocksY,derivsAY,derivsBY,derivsCY,derivnRHSy,derivsRY)
        case(1)
            call solveSystemYper(F2,Fy2,nDerivBlocksY,derivBlocksY,derivsAY,derivsBY,derivsCY,derivsDY,derivnRHSy,derivsRY)
        end select
        call transpose_y_to_x(Fy2,Fy)
            
    end subroutine
    
    subroutine calcDerivsXZ(F,Fx,Fz,nDerivBlocksX,nDerivBlocksZ,derivBlocksX,derivBlocksZ, &
                        & derivnRHSx,derivnRHSz, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDX, derivsDZ, &
                        & periodicX, periodicZ)
                        
        ! INPUTS
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: F
        integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksZ
        integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksZ
        integer, intent(in)                   :: derivnRHSx, derivnRHSz
        real*8,  dimension(:,:), intent(in)   :: derivsAX, derivsBX, derivsCX, derivsAZ, derivsBZ, derivsCZ
        real*8,  dimension(:,:), intent(in)   :: derivsDX, derivsDZ
        real*8,  dimension(:,:,:), intent(in) :: derivsRX, derivsRZ
        integer, intent(in)                   :: periodicX, periodicZ
        
        ! OUTPUT
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Fx, Fz
        
        ! INTERNALS
        real*8,  dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: F2, Fz2
        real*8,  dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: F3, Fz3
        
        select case(periodicX)
        case(0)
            call solveSystemX(F,Fx,nDerivBlocksX,derivBlocksX,derivsAX,derivsBX,derivsCX,derivnRHSx,derivsRX)
        case(1)
            call solveSystemXper(F,Fx,nDerivBlocksX,derivBlocksX,derivsAX,derivsBX,derivsCX,derivsDX,derivnRHSx,derivsRX)
        end select
        
        call transpose_x_to_y(F,F2)
        call transpose_y_to_z(F2,F3)
        
        select case(periodicZ)
        case(0)
            call solveSystemZ(F3,Fz3,nDerivBlocksZ,derivBlocksZ,derivsAZ,derivsBZ,derivsCZ,derivnRHSz,derivsRZ)
        case(1)
            call solveSystemZper(F3,Fz3,nDerivBlocksZ,derivBlocksZ,derivsAZ,derivsBZ,derivsCZ,derivsDZ,derivnRHSz,derivsRZ)
        end select
        
        call transpose_z_to_y(Fz3,Fz2)
        call transpose_y_to_x(Fz2,Fz)
            
    end subroutine
    
    subroutine calcDerivsYZ(F,Fy,Fz,nDerivBlocksY,nDerivBlocksZ,derivBlocksY,derivBlocksZ, &
                        & derivnRHSy,derivnRHSz, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDY, derivsDZ, &
                        & periodicY, periodicZ)
                        
        ! INPUTS
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: F
        integer, intent(in)                   :: nDerivBlocksY,nDerivBlocksZ
        integer, dimension(:,:), intent(in)   :: derivBlocksY,derivBlocksZ
        integer, intent(in)                   :: derivnRHSy, derivnRHSz
        real*8,  dimension(:,:), intent(in)   :: derivsAY, derivsBY, derivsCY, derivsAZ, derivsBZ, derivsCZ
        real*8,  dimension(:,:), intent(in)   :: derivsDY, derivsDZ
        real*8,  dimension(:,:,:), intent(in) :: derivsRY, derivsRZ
        integer, intent(in)                   :: periodicY, periodicZ
        
        ! OUTPUT
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Fy, Fz
        
        ! INTERNALS
        real*8,  dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: F2, Fy2, Fz2
        real*8,  dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: F3, Fz3
        
        call transpose_x_to_y(F,F2)
        select case(periodicY)
        case(0)
            call solveSystemY(F2,Fy2,nDerivBlocksY,derivBlocksY,derivsAY,derivsBY,derivsCY,derivnRHSy,derivsRY)
        case(1)
            call solveSystemYper(F2,Fy2,nDerivBlocksY,derivBlocksY,derivsAY,derivsBY,derivsCY,derivsDY,derivnRHSy,derivsRY)
        end select
        call transpose_y_to_x(Fy2,Fy)
        
        call transpose_y_to_z(F2,F3)
        
        select case(periodicZ)
        case(0)
            call solveSystemZ(F3,Fz3,nDerivBlocksZ,derivBlocksZ,derivsAZ,derivsBZ,derivsCZ,derivnRHSz,derivsRZ)
        case(1)
            call solveSystemZper(F3,Fz3,nDerivBlocksZ,derivBlocksZ,derivsAZ,derivsBZ,derivsCZ,derivsDZ,derivnRHSz,derivsRZ)
        end select
        
        call transpose_z_to_y(Fz3,Fz2)
        call transpose_y_to_x(Fz2,Fz)
            
    end subroutine
    
    subroutine calcDerivsXYZ(F,Fx,Fy,Fz,nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                        & derivnRHSx,derivnRHSy,derivnRHSz, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDX, derivsDY, derivsDZ, &
                        & periodicX, periodicY, periodicZ)
                        
        ! INPUTS
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: F
        integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ
        integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksY,derivBlocksZ
        integer, intent(in)                   :: derivnRHSx, derivnRHSy, derivnRHSz
        real*8,  dimension(:,:), intent(in)   :: derivsAX, derivsBX, derivsCX, derivsAY, derivsBY, derivsCY, derivsAZ, derivsBZ, derivsCZ
        real*8,  dimension(:,:), intent(in)   :: derivsDX, derivsDY, derivsDZ
        real*8,  dimension(:,:,:), intent(in) :: derivsRX, derivsRY, derivsRZ
        integer, intent(in)                   :: periodicX, periodicY, periodicZ
        
        ! OUTPUT
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Fx, Fy, Fz
        
        ! INTERNALS
        real*8,  dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: F2, Fy2, Fz2
        real*8,  dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: F3, Fz3
        
        select case(periodicX)
        case(0)
            call solveSystemX(F,Fx,nDerivBlocksX,derivBlocksX,derivsAX,derivsBX,derivsCX,derivnRHSx,derivsRX)
        case(1)
            call solveSystemXper(F,Fx,nDerivBlocksX,derivBlocksX,derivsAX,derivsBX,derivsCX,derivsDX,derivnRHSx,derivsRX)
        end select
        
        call transpose_x_to_y(F,F2)
        select case(periodicY)
        case(0)
            call solveSystemY(F2,Fy2,nDerivBlocksY,derivBlocksY,derivsAY,derivsBY,derivsCY,derivnRHSy,derivsRY)
        case(1)
            call solveSystemYper(F2,Fy2,nDerivBlocksY,derivBlocksY,derivsAY,derivsBY,derivsCY,derivsDY,derivnRHSy,derivsRY)
        end select
        call transpose_y_to_x(Fy2,Fy)
        
        call transpose_y_to_z(F2,F3)
        
        select case(periodicZ)
        case(0)
            call solveSystemZ(F3,Fz3,nDerivBlocksZ,derivBlocksZ,derivsAZ,derivsBZ,derivsCZ,derivnRHSz,derivsRZ)
        case(1)
            call solveSystemZper(F3,Fz3,nDerivBlocksZ,derivBlocksZ,derivsAZ,derivsBZ,derivsCZ,derivsDZ,derivnRHSz,derivsRZ)
        end select
        
        call transpose_z_to_y(Fz3,Fz2)
        call transpose_y_to_x(Fz2,Fz)
            
    end subroutine
    
    ! APPLY FILTERS
    
    subroutine applyFilters2D(U,V,R,E,FilterX,FilterY, &
                        & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                        & filternRHSx,filternRHSy, &
                        & filterAX, filterBX, filterCX, filterRX, &
                        & filterAY, filterBY, filterCY, filterRY, &
                        & filterDX, filterDY, &
                        & periodicX, periodicY, &
						& FilterCharTime)
                        
                        
        ! INPUTS
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(inout) :: U,V,R,E
        integer, intent(in)                   :: FilterX, FilterY
        integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksY
        integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksY
        integer, intent(in)                   :: filternRHSx, filternRHSy
        real*8,  dimension(:,:), intent(in)   :: filterAX, filterBX, filterCX, filterAY, filterBY, filterCY
        real*8,  dimension(:,:), intent(in)   :: filterDX, filterDY
        real*8,  dimension(:,:,:), intent(in) :: filterRX, filterRY
        integer, intent(in)                   :: periodicX, periodicY
		real*8, intent(in)                    :: FilterCharTime
        
        ! INTERNALS
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: Uf, Vf, Rf, Ef
        real*8,  dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: U2, V2, R2, E2, Uf2, Vf2, Rf2, Ef2

        
        if(FilterX.eq.1) then
            select case(periodicX)
            case(0)
                call solveSystemX(U,Uf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filternRHSx,filterRX)
                call solveSystemX(V,Vf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filternRHSx,filterRX)
                call solveSystemX(R,Rf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filternRHSx,filterRX)
                call solveSystemX(E,Ef,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filternRHSx,filterRX)
            case(1)
                call solveSystemXper(U,Uf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filterDX,filternRHSx,filterRX)
                call solveSystemXper(V,Vf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filterDX,filternRHSx,filterRX)
                call solveSystemXper(R,Rf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filterDX,filternRHSx,filterRX)
                call solveSystemXper(E,Ef,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filterDX,filternRHSx,filterRX)
            end select
			
			if(FilterCharTime.gt.0) then
				Uf = U + FilterCharTime * (Uf-U)
				Vf = V + FilterCharTime * (Vf-V)
				Rf = R + FilterCharTime * (Rf-R)
				Ef = E + FilterCharTime * (Ef-E)
			endif
			
        else
            Uf = U
            Vf = V
            Rf = R
            Ef = E
        endif
        
        if(FilterY.eq.1) then
        
            call transpose_x_to_y(Uf,U2)
            call transpose_x_to_y(Vf,V2)
            call transpose_x_to_y(Rf,R2)
            call transpose_x_to_y(Ef,E2)
        
            select case(periodicY)
            case(0)
                call solveSystemY(U2,Uf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filternRHSy,filterRY)
                call solveSystemY(V2,Vf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filternRHSy,filterRY)
                call solveSystemY(R2,Rf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filternRHSy,filterRY)
                call solveSystemY(E2,Ef2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filternRHSy,filterRY)
            case(1)
                call solveSystemYper(U2,Uf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filterDY,filternRHSy,filterRY)
                call solveSystemYper(V2,Vf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filterDY,filternRHSy,filterRY)
                call solveSystemYper(R2,Rf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filterDY,filternRHSy,filterRY)
                call solveSystemYper(E2,Ef2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filterDY,filternRHSy,filterRY)
            end select
			
			if(FilterCharTime.gt.0) then
				Uf2 = U2 + FilterCharTime * (Uf2-U2)
				Vf2 = V2 + FilterCharTime * (Vf2-V2)
				Rf2 = R2 + FilterCharTime * (Rf2-R2)
				Ef2 = E2 + FilterCharTime * (Ef2-E2)
			endif
            
            call transpose_y_to_x(Uf2,U)
            call transpose_y_to_x(Vf2,V)
            call transpose_y_to_x(Rf2,R)
            call transpose_y_to_x(Ef2,E)
            
        else
            U = Uf
            V = Vf
            R = Rf
            E = Ef
        endif
        
    end subroutine

    subroutine applyFilters3D(U,V,W,R,E,FilterX,FilterY,FilterZ, &
                        & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                        & filternRHSx,filternRHSy,filternRHSz, &
                        & filterAX, filterBX, filterCX, filterRX, &
                        & filterAY, filterBY, filterCY, filterRY, &
                        & filterAZ, filterBZ, filterCZ, filterRZ, &
                        & filterDX, filterDY, filterDZ, &
                        & periodicX, periodicY, periodicZ, &
						& FilterCharTime)
                        
                        
        ! INPUTS
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(inout) :: U,V,W,R,E
        integer, intent(in)                   :: FilterX, FilterY, FilterZ
        integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksY, nDerivBlocksZ
        integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksY, derivBlocksZ
        integer, intent(in)                   :: filternRHSx, filternRHSy, filternRHSz
        real*8,  dimension(:,:), intent(in)   :: filterAX, filterBX, filterCX, filterAY, filterBY, filterCY, filterAZ, filterBZ, filterCZ
        real*8,  dimension(:,:), intent(in)   :: filterDX, filterDY, filterDZ
        real*8,  dimension(:,:,:), intent(in) :: filterRX, filterRY, filterRZ
        integer, intent(in)                   :: periodicX, periodicY, periodicZ
		real*8, intent(in)                    :: FilterCharTime
        
        ! INTERNALS
        real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: Uf, Vf, Wf, Rf, Ef
        real*8,  dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: U2, V2, W2, R2, E2, Uf2, Vf2, Wf2, Rf2, Ef2
        real*8,  dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: U3, V3, W3, R3, E3, Uf3, Vf3, Wf3, Rf3, Ef3
        
        if(FilterX.eq.1) then
            select case(periodicX)
            case(0)
                call solveSystemX(U,Uf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filternRHSx,filterRX)
                call solveSystemX(V,Vf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filternRHSx,filterRX)
                call solveSystemX(W,Wf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filternRHSx,filterRX)
                call solveSystemX(R,Rf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filternRHSx,filterRX)
                call solveSystemX(E,Ef,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filternRHSx,filterRX)
            case(1)
                call solveSystemXper(U,Uf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filterDX,filternRHSx,filterRX)
                call solveSystemXper(V,Vf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filterDX,filternRHSx,filterRX)
                call solveSystemXper(W,Wf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filterDX,filternRHSx,filterRX)
                call solveSystemXper(R,Rf,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filterDX,filternRHSx,filterRX)
                call solveSystemXper(E,Ef,nDerivBlocksX,derivBlocksX,filterAX,filterBX,filterCX,filterDX,filternRHSx,filterRX)
            end select
			
			if(FilterCharTime.gt.0) then
				Uf = U + FilterCharTime * (Uf-U)
				Vf = V + FilterCharTime * (Vf-V)
				Wf = W + FilterCharTime * (Wf-W)
				Rf = R + FilterCharTime * (Rf-R)
				Ef = E + FilterCharTime * (Ef-E)
			endif
			
        else
            Uf = U
            Vf = V
            Wf = W
            Rf = R
            Ef = E
        endif
        
        if((FilterY.eq.1).or.(FilterZ.eq.1)) then
        
            call transpose_x_to_y(Uf,U2)
            call transpose_x_to_y(Vf,V2)
            call transpose_x_to_y(Wf,W2)
            call transpose_x_to_y(Rf,R2)
            call transpose_x_to_y(Ef,E2)
        
            if(FilterY.eq.1) then
                select case(periodicY)
                case(0)
                    call solveSystemY(U2,Uf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filternRHSy,filterRY)
                    call solveSystemY(V2,Vf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filternRHSy,filterRY)
                    call solveSystemY(W2,Wf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filternRHSy,filterRY)
                    call solveSystemY(R2,Rf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filternRHSy,filterRY)
                    call solveSystemY(E2,Ef2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filternRHSy,filterRY)
                case(1)
                    call solveSystemYper(U2,Uf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filterDY,filternRHSy,filterRY)
                    call solveSystemYper(V2,Vf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filterDY,filternRHSy,filterRY)
                    call solveSystemYper(W2,Wf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filterDY,filternRHSy,filterRY)
                    call solveSystemYper(R2,Rf2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filterDY,filternRHSy,filterRY)
                    call solveSystemYper(E2,Ef2,nDerivBlocksY,derivBlocksY,filterAY,filterBY,filterCY,filterDY,filternRHSy,filterRY)
                end select
				
				if(FilterCharTime.gt.0) then
					Uf2 = U2 + FilterCharTime * (Uf2-U2)
					Vf2 = V2 + FilterCharTime * (Vf2-V2)
					Wf2 = W2 + FilterCharTime * (Wf2-W2)
					Rf2 = R2 + FilterCharTime * (Rf2-R2)
					Ef2 = E2 + FilterCharTime * (Ef2-E2)
				endif
				
            else
                Uf2 = U2
                Vf2 = V2
                Wf2 = W2
                Rf2 = R2
                Ef2 = E2
            endif
            
            if(FilterZ.eq.1) then
            
                call transpose_y_to_z(Uf2,U3)
                call transpose_y_to_z(Vf2,V3)
                call transpose_y_to_z(Wf2,W3)
                call transpose_y_to_z(Rf2,R3)
                call transpose_y_to_z(Ef2,E3)
            
                select case(periodicZ)
                case(0)
                    call solveSystemZ(U3,Uf3,nDerivBlocksZ,derivBlocksZ,filterAZ,filterBZ,filterCZ,filternRHSz,filterRZ)
                    call solveSystemZ(V3,Vf3,nDerivBlocksZ,derivBlocksZ,filterAZ,filterBZ,filterCZ,filternRHSz,filterRZ)
                    call solveSystemZ(W3,Wf3,nDerivBlocksZ,derivBlocksZ,filterAZ,filterBZ,filterCZ,filternRHSz,filterRZ)
                    call solveSystemZ(R3,Rf3,nDerivBlocksZ,derivBlocksZ,filterAZ,filterBZ,filterCZ,filternRHSz,filterRZ)
                    call solveSystemZ(E3,Ef3,nDerivBlocksZ,derivBlocksZ,filterAZ,filterBZ,filterCZ,filternRHSz,filterRZ)
                case(1)
                    call solveSystemZper(U3,Uf3,nDerivBlocksZ,derivBlocksZ,filterAZ,filterBZ,filterCZ,filterDZ,filternRHSz,filterRZ)
                    call solveSystemZper(V3,Vf3,nDerivBlocksZ,derivBlocksZ,filterAZ,filterBZ,filterCZ,filterDZ,filternRHSz,filterRZ)
                    call solveSystemZper(W3,Wf3,nDerivBlocksZ,derivBlocksZ,filterAZ,filterBZ,filterCZ,filterDZ,filternRHSz,filterRZ)
                    call solveSystemZper(R3,Rf3,nDerivBlocksZ,derivBlocksZ,filterAZ,filterBZ,filterCZ,filterDZ,filternRHSz,filterRZ)
                    call solveSystemZper(E3,Ef3,nDerivBlocksZ,derivBlocksZ,filterAZ,filterBZ,filterCZ,filterDZ,filternRHSz,filterRZ)
                end select
				
				if(FilterCharTime.gt.0) then
					Uf3 = U3 + FilterCharTime * (Uf3-U3)
					Vf3 = V3 + FilterCharTime * (Vf3-V3)
					Wf3 = W3 + FilterCharTime * (Wf3-W3)
					Rf3 = R3 + FilterCharTime * (Rf3-R3)
					Ef3 = E3 + FilterCharTime * (Ef3-E3)
				endif
                
                call transpose_z_to_y(Uf3,Uf2)
                call transpose_z_to_y(Vf3,Vf2)
                call transpose_z_to_y(Wf3,Wf2)
                call transpose_z_to_y(Rf3,Rf2)
                call transpose_z_to_y(Ef3,Ef2)
                
            endif
            
            call transpose_y_to_x(Uf2,U)
            call transpose_y_to_x(Vf2,V)
            call transpose_y_to_x(Wf2,W)
            call transpose_y_to_x(Rf2,R)
            call transpose_y_to_x(Ef2,E)
            
        else
            U = Uf
            V = Vf
            W = Wf
            R = Rf
            E = Ef
        endif
        
    end subroutine
    
    ! SOLVE THE SYSTEMS
    
    subroutine solveSystemX(F,dF,nBlocks,blocks,A,B,C,nR,R)
    
    ! INPUTS
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: F
    integer, intent(in)                   :: nBlocks, nR
    integer, dimension(:,:), intent(in)   :: blocks
    real*8,  dimension(:,:), intent(in)   :: A,B,C
    real*8,  dimension(:,:,:), intent(in) :: R
    
    ! OUTPUT
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: dF
    
    ! INTERNALS
    integer :: i,l,bn,t
    integer :: j1,j2,k1,k2
    
    do bn = 1,nBlocks
    
        ! Get block indices
        t = blocks(bn,1)
        
        j1 = blocks(bn,2)
        j2 = blocks(bn,3)
        k1 = blocks(bn,4)
        k2 = blocks(bn,5)
        
        ! Compute the right hand side
        do i = 1,nR
            dF(i,j1:j2,k1:k2) = R(i,nR,t) * F(i,j1:j2,k1:k2)
            do l = 1,nR-1
                dF(i,j1:j2,k1:k2) = dF(i,j1:j2,k1:k2) + R(i,nR+l,t) * F(i+l,j1:j2,k1:k2) + R(i,nR-l,t) * F(1+modulo(i-1-l,xsize(1)),j1:j2,k1:k2)
            end do
        end do
        
        do i = nR+1,xsize(1)-nR
            dF(i,j1:j2,k1:k2) = R(i,nR,t) * F(i,j1:j2,k1:k2)
            do l = 1,nR-1
                dF(i,j1:j2,k1:k2) = dF(i,j1:j2,k1:k2) + R(i,nR+l,t) * F(i+l,j1:j2,k1:k2) + R(i,nR-l,t) * F(i-l,j1:j2,k1:k2)
            end do
        end do
        
        do i = xsize(1)-nR+1,xsize(1)
            dF(i,j1:j2,k1:k2) = R(i,nR,t) * F(i,j1:j2,k1:k2)
            do l = 1,nR-1
                dF(i,j1:j2,k1:k2) = dF(i,j1:j2,k1:k2) + R(i,nR+l,t) * F(1+modulo(i-1+l,xsize(1)),j1:j2,k1:k2) + R(i,nR-l,t) * F(i-l,j1:j2,k1:k2)
            end do
        end do
        
        ! Solve the left hand side
        dF(1,j1:j2,k1:k2) = dF(1,j1:j2,k1:k2)*B(1,t)
        do i = 1,xsize(1)-1
            dF(i+1,j1:j2,k1:k2) = (dF(i+1,j1:j2,k1:k2) - A(i,t)*dF(i,j1:j2,k1:k2))*B(i+1,t)
        end do
        do i = xsize(1)-1,1,-1
            dF(i,j1:j2,k1:k2) = dF(i,j1:j2,k1:k2) - C(i,t)*dF(i+1,j1:j2,k1:k2)
        end do
        
    end do
    
    end subroutine
    
    subroutine solveSystemY(F,dF,nBlocks,blocks,A,B,C,nR,R)
    
    ! INPUTS
    real*8,  dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), intent(in) :: F
    integer, intent(in)                   :: nBlocks, nR
    integer, dimension(:,:), intent(in)   :: blocks
    real*8,  dimension(:,:), intent(in)   :: A,B,C
    real*8,  dimension(:,:,:), intent(in) :: R
    
    ! OUTPUT
    real*8,  dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), intent(out) :: dF
    
    ! INTERNALS
    integer :: j,l,bn,t
    integer :: i1,i2,k1,k2
    
    do bn = 1,nBlocks
    
        ! Get block indices
        t = blocks(bn,1)
        
        i1 = blocks(bn,2)
        i2 = blocks(bn,3)
        k1 = blocks(bn,4)
        k2 = blocks(bn,5)
        
        ! Compute the right hand side
        do j = 1,nR
            dF(i1:i2,j,k1:k2) = R(j,nR,t) * F(i1:i2,j,k1:k2)
            do l = 1,nR-1
                dF(i1:i2,j,k1:k2) = dF(i1:i2,j,k1:k2) + R(j,nR+l,t) * F(i1:i2,j+l,k1:k2) + R(j,nR-l,t) * F(i1:i2,1+modulo(j-1-l,ysize(2)),k1:k2)
            end do
        end do
        
        do j = nR+1,ysize(2)-nR
            dF(i1:i2,j,k1:k2) = R(j,nR,t) * F(i1:i2,j,k1:k2)
            do l = 1,nR-1
                dF(i1:i2,j,k1:k2) = dF(i1:i2,j,k1:k2) + R(j,nR+l,t) * F(i1:i2,j+l,k1:k2) + R(j,nR-l,t) * F(i1:i2,j-l,k1:k2)
            end do
        end do
        
        do j = ysize(2)-nR+1,ysize(2)
            dF(i1:i2,j,k1:k2) = R(j,nR,t) * F(i1:i2,j,k1:k2)
            do l = 1,nR-1
                dF(i1:i2,j,k1:k2) = dF(i1:i2,j,k1:k2) + R(j,nR+l,t) * F(i1:i2,1+modulo(j-1+l,ysize(2)),k1:k2) + R(j,nR-l,t) * F(i1:i2,j-l,k1:k2)
            end do
        end do
        
        ! Solve the left hand side
        dF(i1:i2,1,k1:k2) = dF(i1:i2,1,k1:k2)*B(1,t)
        do j = 1,ysize(2)-1
            dF(i1:i2,j+1,k1:k2) = (dF(i1:i2,j+1,k1:k2) - A(j,t)*dF(i1:i2,j,k1:k2))*B(j+1,t)
        end do
        do j = ysize(2)-1,1,-1
            dF(i1:i2,j,k1:k2) = dF(i1:i2,j,k1:k2) - C(j,t)*dF(i1:i2,j+1,k1:k2)
        end do
    end do
    
    end subroutine
    
    subroutine solveSystemZ(F,dF,nBlocks,blocks,A,B,C,nR,R)
    
    ! INPUTS
    real*8,  dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)), intent(in) :: F
    integer, intent(in)                   :: nBlocks, nR
    integer, dimension(:,:), intent(in)   :: blocks
    real*8,  dimension(:,:), intent(in)   :: A,B,C
    real*8,  dimension(:,:,:), intent(in) :: R
    
    ! OUTPUT
    real*8,  dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)), intent(out) :: dF
    
    ! INTERNALS
    integer :: k,l,bn,t
    integer :: i1,i2,j1,j2
    
    do bn = 1,nBlocks
    
        ! Get block indices
        t = blocks(bn,1)
        
        i1 = blocks(bn,2)
        i2 = blocks(bn,3)
        j1 = blocks(bn,4)
        j2 = blocks(bn,5)
        
        ! Compute the right hand side
        do k = 1,nR
            dF(i1:i2,j1:j2,k) = R(k,nR,t) * F(i1:i2,j1:j2,k)
            do l = 1,nR-1
                dF(i1:i2,j1:j2,k) = dF(i1:i2,j1:j2,k) + R(k,nR+l,t) * F(i1:i2,j1:j2,k+l) + R(k,nR-l,t) * F(i1:i2,j1:j2,1+modulo(k-1-l,zsize(3)))
            end do
        end do
        
        do k = nR+1,zsize(3)-nR
            dF(i1:i2,j1:j2,k) = R(k,nR,t) * F(i1:i2,j1:j2,k)
            do l = 1,nR-1
                dF(i1:i2,j1:j2,k) = dF(i1:i2,j1:j2,k) + R(k,nR+l,t) * F(i1:i2,j1:j2,k+l) + R(k,nR-l,t) * F(i1:i2,j1:j2,k-l)
            end do
        end do
        
        do k = zsize(3)-nR+1,zsize(3)
            dF(i1:i2,j1:j2,k) = R(k,nR,t) * F(i1:i2,j1:j2,k)
            do l = 1,nR-1
                dF(i1:i2,j1:j2,k) = dF(i1:i2,j1:j2,k) + R(k,nR+l,t) * F(i1:i2,j1:j2,1+modulo(k-1+l,zsize(3))) + R(k,nR-l,t) * F(i1:i2,j1:j2,k-l)
            end do
        end do
        
        ! Solve the left hand side
        dF(i1:i2,j1:j2,1) = dF(i1:i2,j1:j2,1)*B(1,t)
        do k = 1,zsize(3)-1
            dF(i1:i2,j1:j2,k+1) = (dF(i1:i2,j1:j2,k+1) - A(k,t)*dF(i1:i2,j1:j2,k))*B(k+1,t)
        end do
        do k = zsize(3)-1,1,-1
            dF(i1:i2,j1:j2,k) = dF(i1:i2,j1:j2,k) - C(k,t)*dF(i1:i2,j1:j2,k+1)
        end do
    end do
    
    end subroutine
    
    subroutine solveSystemXper(F,dF,nBlocks,blocks,A,B,C,D,nR,R)
    
    ! INPUTS
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: F
    integer, intent(in)                   :: nBlocks, nR
    integer, dimension(:,:), intent(in)   :: blocks
    real*8,  dimension(:,:), intent(in)   :: A,B,C,D
    real*8,  dimension(:,:,:), intent(in) :: R
    
    ! OUTPUT
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: dF
    
    ! INTERNALS
    integer :: i,l,bn,t
    integer :: j1,j2,k1,k2
    
    do bn = 1,nBlocks
    
        ! Get block indices
        t = blocks(bn,1)
        
        j1 = blocks(bn,2)
        j2 = blocks(bn,3)
        k1 = blocks(bn,4)
        k2 = blocks(bn,5)
        
        ! Compute the right hand side
        do i = 1,nR
            dF(i,j1:j2,k1:k2) = R(i,nR,t) * F(i,j1:j2,k1:k2)
            do l = 1,nR-1
                dF(i,j1:j2,k1:k2) = dF(i,j1:j2,k1:k2) + R(i,nR+l,t) * F(i+l,j1:j2,k1:k2) + R(i,nR-l,t) * F(1+modulo(i-1-l,xsize(1)),j1:j2,k1:k2)
            end do
        end do
        
        do i = nR+1,xsize(1)-nR
            dF(i,j1:j2,k1:k2) = R(i,nR,t) * F(i,j1:j2,k1:k2)
            do l = 1,nR-1
                dF(i,j1:j2,k1:k2) = dF(i,j1:j2,k1:k2) + R(i,nR+l,t) * F(i+l,j1:j2,k1:k2) + R(i,nR-l,t) * F(i-l,j1:j2,k1:k2)
            end do
        end do
        
        do i = xsize(1)-nR+1,xsize(1)
            dF(i,j1:j2,k1:k2) = R(i,nR,t) * F(i,j1:j2,k1:k2)
            do l = 1,nR-1
                dF(i,j1:j2,k1:k2) = dF(i,j1:j2,k1:k2) + R(i,nR+l,t) * F(1+modulo(i-1+l,xsize(1)),j1:j2,k1:k2) + R(i,nR-l,t) * F(i-l,j1:j2,k1:k2)
            end do
        end do
        
        ! Solve the first element of the left hand side
        dF(1,j1:j2,k1:k2) = dF(1,j1:j2,k1:k2)*D(1,t)
        do i = 2,xsize(1)
            dF(1,j1:j2,k1:k2) = dF(1,j1:j2,k1:k2) + dF(i,j1:j2,k1:k2)*D(i,t)
        end do
        
        ! Add the known solution to the right hand side
        dF(2,j1:j2,k1:k2) = dF(2,j1:j2,k1:k2) - A(1,t) * dF(1,j1:j2,k1:k2)
        dF(xsize(1),j1:j2,k1:k2) = dF(xsize(1),j1:j2,k1:k2) - C(1,t) * dF(1,j1:j2,k1:k2)
        
        ! Solve the left hand side
        dF(2,j1:j2,k1:k2) = dF(2,j1:j2,k1:k2)*B(1,t)
        do i = 2,xsize(1)-1
            dF(i+1,j1:j2,k1:k2) = (dF(i+1,j1:j2,k1:k2) - A(i,t)*dF(i,j1:j2,k1:k2))*B(i,t)
        end do
        do i = xsize(1)-1,2,-1
            dF(i,j1:j2,k1:k2) = dF(i,j1:j2,k1:k2) - C(i,t)*dF(i+1,j1:j2,k1:k2)
        end do
        
    end do
    
    end subroutine
    
    subroutine solveSystemYper(F,dF,nBlocks,blocks,A,B,C,D,nR,R)
    
    ! INPUTS
    real*8,  dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), intent(in) :: F
    integer, intent(in)                   :: nBlocks, nR
    integer, dimension(:,:), intent(in)   :: blocks
    real*8,  dimension(:,:), intent(in)   :: A,B,C,D
    real*8,  dimension(:,:,:), intent(in) :: R
    
    ! OUTPUT
    real*8,  dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), intent(out) :: dF
    
    ! INTERNALS
    integer :: j,l,bn,t
    integer :: i1,i2,k1,k2
    
    do bn = 1,nBlocks
    
        ! Get block indices
        t = blocks(bn,1)
        
        i1 = blocks(bn,2)
        i2 = blocks(bn,3)
        k1 = blocks(bn,4)
        k2 = blocks(bn,5)
        
        ! Compute the right hand side
        do j = 1,nR
            dF(i1:i2,j,k1:k2) = R(j,nR,t) * F(i1:i2,j,k1:k2)
            do l = 1,nR-1
                dF(i1:i2,j,k1:k2) = dF(i1:i2,j,k1:k2) + R(j,nR+l,t) * F(i1:i2,j+l,k1:k2) + R(j,nR-l,t) * F(i1:i2,1+modulo(j-1-l,ysize(2)),k1:k2)
            end do
        end do
        
        do j = nR+1,ysize(2)-nR
            dF(i1:i2,j,k1:k2) = R(j,nR,t) * F(i1:i2,j,k1:k2)
            do l = 1,nR-1
                dF(i1:i2,j,k1:k2) = dF(i1:i2,j,k1:k2) + R(j,nR+l,t) * F(i1:i2,j+l,k1:k2) + R(j,nR-l,t) * F(i1:i2,j-l,k1:k2)
            end do
        end do
        
        do j = ysize(2)-nR+1,ysize(2)
            dF(i1:i2,j,k1:k2) = R(j,nR,t) * F(i1:i2,j,k1:k2)
            do l = 1,nR-1
                dF(i1:i2,j,k1:k2) = dF(i1:i2,j,k1:k2) + R(j,nR+l,t) * F(i1:i2,1+modulo(j-1+l,ysize(2)),k1:k2) + R(j,nR-l,t) * F(i1:i2,j-l,k1:k2)
            end do
        end do
        
        ! Solve the first element of the left hand side
        dF(i1:i2,1,k1:k2) = dF(i1:i2,1,k1:k2)*D(1,t)
        do j = 2,ysize(2)
            dF(i1:i2,1,k1:k2) = dF(i1:i2,1,k1:k2) + dF(i1:i2,j,k1:k2)*D(j,t)
        end do
        
        ! Add the known solution to the right hand side
        dF(i1:i2,2,k1:k2) = dF(i1:i2,2,k1:k2) - A(1,t) * dF(i1:i2,1,k1:k2)
        dF(i1:i2,ysize(2),k1:k2) = dF(i1:i2,ysize(2),k1:k2) - C(1,t) * dF(i1:i2,1,k1:k2)
        
        ! Solve the left hand side
        dF(i1:i2,2,k1:k2) = dF(i1:i2,2,k1:k2)*B(1,t)
        do j = 2,ysize(2)-1
            dF(i1:i2,j+1,k1:k2) = (dF(i1:i2,j+1,k1:k2) - A(j,t)*dF(i1:i2,j,k1:k2))*B(j,t)
        end do
        do j = ysize(2)-1,2,-1
            dF(i1:i2,j,k1:k2) = dF(i1:i2,j,k1:k2) - C(j,t)*dF(i1:i2,j+1,k1:k2)
        end do
        
    end do
    
    end subroutine
    
    subroutine solveSystemZper(F,dF,nBlocks,blocks,A,B,C,D,nR,R)
    
    ! INPUTS
    real*8,  dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)), intent(in) :: F
    integer, intent(in)                   :: nBlocks, nR
    integer, dimension(:,:), intent(in)   :: blocks
    real*8,  dimension(:,:), intent(in)   :: A,B,C,D
    real*8,  dimension(:,:,:), intent(in) :: R
    
    ! OUTPUT
    real*8,  dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)), intent(out) :: dF
    
    ! INTERNALS
    integer :: k,l,bn,t
    integer :: i1,i2,j1,j2
    
    do bn = 1,nBlocks
    
        ! Get block indices
        t = blocks(bn,1)
        
        i1 = blocks(bn,2)
        i2 = blocks(bn,3)
        j1 = blocks(bn,4)
        j2 = blocks(bn,5)
        
        ! Compute the right hand side
        do k = 1,nR
            dF(i1:i2,j1:j2,k) = R(k,nR,t) * F(i1:i2,j1:j2,k)
            do l = 1,nR-1
                dF(i1:i2,j1:j2,k) = dF(i1:i2,j1:j2,k) + R(k,nR+l,t) * F(i1:i2,j1:j2,k+l) + R(k,nR-l,t) * F(i1:i2,j1:j2,1+modulo(k-1-l,zsize(3)))
            end do
        end do
        
        do k = nR+1,zsize(3)-nR
            dF(i1:i2,j1:j2,k) = R(k,nR,t) * F(i1:i2,j1:j2,k)
            do l = 1,nR-1
                dF(i1:i2,j1:j2,k) = dF(i1:i2,j1:j2,k) + R(k,nR+l,t) * F(i1:i2,j1:j2,k+l) + R(k,nR-l,t) * F(i1:i2,j1:j2,k-l)
            end do
        end do
        
        do k = zsize(3)-nR+1,zsize(3)
            dF(i1:i2,j1:j2,k) = R(k,nR,t) * F(i1:i2,j1:j2,k)
            do l = 1,nR-1
                dF(i1:i2,j1:j2,k) = dF(i1:i2,j1:j2,k) + R(k,nR+l,t) * F(i1:i2,j1:j2,1+modulo(k-1+l,zsize(3))) + R(k,nR-l,t) * F(i1:i2,j1:j2,k-l)
            end do
        end do
        
        ! Solve the first element of the left hand side
        dF(i1:i2,j1:j2,1) = dF(i1:i2,j1:j2,1)*D(1,t)
        do k = 2,zsize(3)
            dF(i1:i2,j1:j2,1) = dF(i1:i2,j1:j2,1) + dF(i1:i2,j1:j2,k)*D(k,t)
        end do
        
        ! Add the known solution to the right hand side
        dF(i1:i2,j1:j2,2) = dF(i1:i2,j1:j2,2) - A(1,t) * dF(i1:i2,j1:j2,1)
        dF(i1:i2,j1:j2,zsize(3)) = dF(i1:i2,j1:j2,zsize(3)) - C(1,t) * dF(i1:i2,j1:j2,1)
        
        ! Solve the left hand side
        dF(i1:i2,j1:j2,2) = dF(i1:i2,j1:j2,2)*B(1,t)
        do k = 2,zsize(3)-1
            dF(i1:i2,j1:j2,k+1) = (dF(i1:i2,j1:j2,k+1) - A(k,t)*dF(i1:i2,j1:j2,k))*B(k,t)
        end do
        do k = zsize(3)-1,2,-1
            dF(i1:i2,j1:j2,k) = dF(i1:i2,j1:j2,k) - C(k,t)*dF(i1:i2,j1:j2,k+1)
        end do
        
    end do
    
    end subroutine
    
    end module
