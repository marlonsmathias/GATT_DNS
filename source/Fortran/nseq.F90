    module nseq
    
    use decomp_2d
    use derivs
    use disturbances
	use mpi
    
    contains
    
    subroutine nseq2d(U,V,R,E,dU,dV,dR,dE, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                & derivnRHSx,derivnRHSy, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsDX, derivsDY, &
                & periodicX, periodicY)
                
    implicit none
    
    ! INPUTS
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: U,V,R,E
    real*8,  intent(in)                   :: Re, Ma, Pr, gamma, T0
    integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksY
    integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksY
    integer, intent(in)                   :: derivnRHSx, derivnRHSy
    real*8,  dimension(:,:), intent(in)   :: derivsAX, derivsBX, derivsCX, derivsAY, derivsBY, derivsCY
    real*8,  dimension(:,:), intent(in)   :: derivsDX, derivsDY
    real*8,  dimension(:,:,:), intent(in) :: derivsRX, derivsRY
    integer, intent(in)                   :: periodicX, periodicY
    
    ! OUTPUTS
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: dU,dV,dR,dE
    
    ! INTERNAL VARIABLES
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: DUx, DVx, DRx, DEx, DUy, DVy, DRy, DEy
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: P, DPx, DPy
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: mu, mu1, mu2, tauxx, tauyy, tauxy
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: dtauxxdx, dtauyydy, dtauxydx, dtauxydy, vtx, vty
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: qx, qy, htx, hty
    
    ! Get spatial derivatives
    
    call calcDerivsXY(U,DUx,DUy,nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                        & derivnRHSx,derivnRHSy, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsDX, derivsDY, &
                        & periodicX, periodicY)
    
    call calcDerivsXY(V,DVx,DVy,nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                        & derivnRHSx,derivnRHSy, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsDX, derivsDY, &
                        & periodicX, periodicY)
    
    call calcDerivsXY(R,DRx,DRy,nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                        & derivnRHSx,derivnRHSy, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsDX, derivsDY, &
                        & periodicX, periodicY)
     
    call calcDerivsXY(E,DEx,DEy,nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                        & derivnRHSx,derivnRHSy, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsDX, derivsDY, &
                        & periodicX, periodicY)
    
    P = (gamma-1.0d0) * R * E
    DPx = (gamma-1.0d0) * (R*DEx + DRx*E)
    DPy = (gamma-1.0d0) * (R*DEy + DRy*E)
    
    ! Viscous term
    
    mu1 = (1.0d0+110.0d0/T0)/(gamma*(gamma-1.0d0)*Ma*Ma*E+110.0d0/T0)
    mu2 = gamma*(gamma-1.0d0)*Ma*Ma*E
    
    mu = mu1*mu2*dsqrt(mu2)
    
    tauxx = 1.0d0/Re * mu * (2.0d0*DUx - 2.0d0/3.0d0 * (DUx+DVy))
    tauyy = 1.0d0/Re * mu * (2.0d0*DVy - 2.0d0/3.0d0 * (DUx+DVy))
    tauxy = 1.0d0/Re * mu * (DUy + DVx)
    
    call calcDerivsXY(tauxy,dtauxydx,dtauxydy,nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                        & derivnRHSx,derivnRHSy, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsDX, derivsDY, &
                        & periodicX, periodicY)
                        
    call calcDerivsX(tauxx,dtauxxdx,nDerivBlocksX,derivBlocksX, &
                        & derivnRHSx, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsDX, &
                        & periodicX)
                        
    call calcDerivsY(tauyy,dtauyydy,nDerivBlocksY,derivBlocksY, &
                        & derivnRHSy, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsDY, &
                        & periodicY)
                        
    vtx = dtauxxdx + dtauxydy
    vty = dtauyydy + dtauxydx
    
    ! Heat term
    
    qx = -gamma/Re/Pr * mu * DEx
    qy = -gamma/Re/Pr * mu * DEy
    
    call calcDerivsX(qx,htx,nDerivBlocksX,derivBlocksX, &
                        & derivnRHSx, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsDX, &
                        & periodicX)
                        
    call calcDerivsY(qy,hty,nDerivBlocksY,derivBlocksY, &
                        & derivnRHSy, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsDY, &
                        & periodicY)

    ! Get derivatives
    
    dR = -R*(DUx+DVy) - DRx*U - DRy*V
    
    dU = -(DUx*U + DUy*V) + (vtx - DPx)/R
    
    dV = -(DVx*U + DVy*V) + (vty - DPy)/R
    
    dE = -(DEx*U + DEy*V) + (-P*(DUx+DVy) + tauxx*DUx + tauyy*DVy + tauxy*(DVx+DUy) - htx - hty)/R
    
    ! CALL FORCINGS
    include 'runForcings2D.F90'

    end subroutine
    
    
    subroutine nseq3d(U,V,W,R,E,dU,dV,dW,dR,dE, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                & derivnRHSx,derivnRHSy,derivnRHSz, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                & derivsDX, derivsDY, derivsDZ, &
                & periodicX, periodicY, periodicZ)
                
    implicit none
    
    ! INPUTS
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: U,V,W,R,E
    real*8,  intent(in)                   :: Re, Ma, Pr, gamma, T0
    integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ
    integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksY,derivBlocksZ
    integer, intent(in)                   :: derivnRHSx, derivnRHSy, derivnRHSz
    real*8,  dimension(:,:), intent(in)   :: derivsAX, derivsBX, derivsCX, derivsAY, derivsBY, derivsCY, derivsAZ, derivsBZ, derivsCZ
    real*8,  dimension(:,:), intent(in)   :: derivsDX, derivsDY, derivsDZ
    real*8,  dimension(:,:,:), intent(in) :: derivsRX, derivsRY, derivsRZ
    integer, intent(in)                   :: periodicX, periodicY, periodicZ
    
    ! OUTPUTS
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: dU,dV,dW,dR,dE
    
    ! INTERNAL VARIABLES
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: DUx, DVx, DWx, DRx, DEx, DUy, DVy, DWy, DRy, DEy, DUz, DVz, DWz, DRz, DEz
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: P, DPx, DPy, DPz
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: mu, mu1, mu2, tauxx, tauyy, tauzz, tauxy, tauxz, tauyz
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: dtauxxdx, dtauyydy, dtauzzdz, vtx, vty, vtz
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: dtauxydx, dtauxydy, dtauxzdx, dtauxzdz, dtauyzdy, dtauyzdz
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: qx, qy, qz, htx, hty, htz
    
    ! Get spatial derivatives
    
    call calcDerivsXYZ(U,DUx,DUy,DUz,nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                        & derivnRHSx,derivnRHSy,derivnRHSz, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDX, derivsDY, derivsDZ, &
                        & periodicX, periodicY, periodicZ)
    
    call calcDerivsXYZ(V,DVx,DVy,DVz,nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                        & derivnRHSx,derivnRHSy,derivnRHSz, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDX, derivsDY, derivsDZ, &
                        & periodicX, periodicY, periodicZ)
                        
    call calcDerivsXYZ(W,DWx,DWy,DWz,nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                        & derivnRHSx,derivnRHSy,derivnRHSz, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDX, derivsDY, derivsDZ, &
                        & periodicX, periodicY, periodicZ)
    
    call calcDerivsXYZ(R,DRx,DRy,DRz,nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                        & derivnRHSx,derivnRHSy,derivnRHSz, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDX, derivsDY, derivsDZ, &
                        & periodicX, periodicY, periodicZ)
    
    call calcDerivsXYZ(E,DEx,DEy,DEz,nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                        & derivnRHSx,derivnRHSy,derivnRHSz, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDX, derivsDY, derivsDZ, &
                        & periodicX, periodicY, periodicZ)
    
    P = (gamma-1.0d0) * R * E
    DPx = (gamma-1.0d0) * (R*DEx + DRx*E)
    DPy = (gamma-1.0d0) * (R*DEy + DRy*E)
    DPz = (gamma-1.0d0) * (R*DEz + DRz*E)
    
    ! Viscous term
    
    mu1 = (1.0d0+110.0d0/T0)/(gamma*(gamma-1.0d0)*Ma*Ma*E+110.0d0/T0)
    mu2 = gamma*(gamma-1.0d0)*Ma*Ma*E
    
    mu = mu1*mu2*dsqrt(mu2)
    
	!tauxx = 1.0d0/Re * mu * (2.0d0*DUx - 2.0d0/3.0d0 * (DUx+DVy+DWz))
    !tauyy = 1.0d0/Re * mu * (2.0d0*DVy - 2.0d0/3.0d0 * (DUx+DVy+DWz))
	!tauzz = 1.0d0/Re * mu * (2.0d0*DWz - 2.0d0/3.0d0 * (DUx+DVy+DWz))
    tauzz = -2.d0/3.d0 * (DUx+DVy+DWz)
    tauxx = 1.0d0/Re * mu * (2.0d0*DUx + tauzz)
    tauyy = 1.0d0/Re * mu * (2.0d0*DVy + tauzz)
    tauzz = 1.0d0/Re * mu * (2.0d0*DWz + tauzz)
    tauxy = 1.0d0/Re * mu * (DUy + DVx)
    tauxz = 1.0d0/Re * mu * (DUz + DWx)
    tauyz = 1.0d0/Re * mu * (DVz + DWy)
    
    call calcDerivsXY(tauxy,dtauxydx,dtauxydy,nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                        & derivnRHSx,derivnRHSy, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsDX, derivsDY, &
                        & periodicX, periodicY)
                        
    call calcDerivsXZ(tauxz,dtauxzdx,dtauxzdz,nDerivBlocksX,nDerivBlocksZ,derivBlocksX,derivBlocksZ, &
                        & derivnRHSx,derivnRHSz, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDX, derivsDZ, &
                        & periodicX, periodicZ)
                        
    call calcDerivsYZ(tauyz,dtauyzdy,dtauyzdz,nDerivBlocksY,nDerivBlocksZ,derivBlocksY,derivBlocksZ, &
                        & derivnRHSy,derivnRHSz, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDY, derivsDZ, &
                        & periodicY, periodicZ)
                        
    call calcDerivsX(tauxx,dtauxxdx,nDerivBlocksX,derivBlocksX, &
                        & derivnRHSx, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsDX, &
                        & periodicX)
                        
    call calcDerivsY(tauyy,dtauyydy,nDerivBlocksY,derivBlocksY, &
                        & derivnRHSy, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsDY, &
                        & periodicY)
                        
   call calcDerivsZ(tauzz,dtauzzdz,nDerivBlocksZ,derivBlocksZ, &
                        & derivnRHSz, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDZ, &
                        & periodicZ)
                        
    vtx = dtauxxdx + dtauxydy + dtauxzdz
    vty = dtauyydy + dtauxydx + dtauyzdz
    vtz = dtauzzdz + dtauxzdx + dtauyzdy
    
    ! Heat term
    
    qx = -gamma/Re/Pr * mu * DEx
    qy = -gamma/Re/Pr * mu * DEy
    qz = -gamma/Re/Pr * mu * DEz
    
    call calcDerivsX(qx,htx,nDerivBlocksX,derivBlocksX, &
                        & derivnRHSx, &
                        & derivsAX, derivsBX, derivsCX, derivsRX, &
                        & derivsDX, &
                        & periodicX)
                        
    call calcDerivsY(qy,hty,nDerivBlocksY,derivBlocksY, &
                        & derivnRHSy, &
                        & derivsAY, derivsBY, derivsCY, derivsRY, &
                        & derivsDY, &
                        & periodicY)
                        
    call calcDerivsZ(qz,htz,nDerivBlocksZ,derivBlocksZ, &
                        & derivnRHSz, &
                        & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                        & derivsDZ, &
                        & periodicZ)

    ! Get derivatives
    
    dR = -R*(DUx+DVy+DWz) - DRx*U - DRy*V - DRz*W
    
    dU = -(DUx*U + DUy*V + DUz*W) + (vtx - DPx)/R
    
    dV = -(DVx*U + DVy*V + DVz*W) + (vty - DPy)/R
    
    dW = -(DWx*U + DWy*V + DWz*W) + (vtz - DPz)/R
    
    dE = -(DEx*U + DEy*V + DEz*W) + (-P*(DUx+DVy+DWz) + tauxx*DUx + tauyy*DVy + tauzz*DWz + tauxy*(DVx+DUy) + tauxz*(DWx+DUz) + tauyz*(DWy+DVz) - htx - hty - htz)/R

    ! CALL FORCINGS
    include 'runForcings3D.F90'
    
    end subroutine
    
    ! BOUNDARY CONDITIONS
    
    subroutine applyBoundCond(U,V,W,R,E, &
                & gamma, t, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    implicit none

    ! INPUT AND OUTPUT
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(inout) :: U,V,W,R,E
    
    ! INPUTS
    real*8,  intent(in)                   :: t, gamma
    integer, intent(in)                   :: neumannLength, neumann2Length
    real*8,  dimension(:), intent(in)     :: neumannCoeffs, neumann2Coeffs
    integer, intent(in)                   :: nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs
    integer, dimension(:,:), intent(in)   :: iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs
    real*8,  dimension(:), intent(in)     :: vUd, vVd, vWd, vPd, vEd
    integer, dimension(:), intent(in)     :: dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs
    integer, intent(in)                   :: cN
    integer, dimension(:,:), intent(in)   :: cL, cD
    integer, dimension(:), intent(in)     :: cAdiabatic
    
    ! INTERNAL VARIABLES
    integer :: i
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: P
    integer, dimension(6) :: ind
    
    ! DIRICHLET CONDITIONS
    
    do i = 1,nUd
        U(iUd(i,1):iUd(i,2),iUd(i,3):iUd(i,4),iUd(i,5):iUd(i,6)) = vUd(i)
    enddo
    do i = 1,nVd
        V(iVd(i,1):iVd(i,2),iVd(i,3):iVd(i,4),iVd(i,5):iVd(i,6)) = vVd(i)
    enddo
    do i = 1,nWd
        W(iWd(i,1):iWd(i,2),iWd(i,3):iWd(i,4),iWd(i,5):iWd(i,6)) = vWd(i)
    enddo
    do i = 1,nEd
        E(iEd(i,1):iEd(i,2),iEd(i,3):iEd(i,4),iEd(i,5):iEd(i,6)) = vEd(i)
    enddo
    do i = 1,nPd
        R(iPd(i,1):iPd(i,2),iPd(i,3):iPd(i,4),iPd(i,5):iPd(i,6)) = vPd(i)/(gamma-1.0d0)/E(iPd(i,1):iPd(i,2),iPd(i,3):iPd(i,4),iPd(i,5):iPd(i,6))
    enddo
    
    
    ! NEUMANN CONDITIONS
    do i = 1,nUn
        call getNeumannIndices(iUn(i,:),ind,dUn(i),neumannLength)
        call applyNeumann(U(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),dUn(i),neumannLength,neumannCoeffs)
    enddo
    do i = 1,nVn
        call getNeumannIndices(iVn(i,:),ind,dVn(i),neumannLength)
        call applyNeumann(V(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),dVn(i),neumannLength,neumannCoeffs)
    enddo
    do i = 1,nWn
        call getNeumannIndices(iWn(i,:),ind,dWn(i),neumannLength)
        call applyNeumann(W(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),dWn(i),neumannLength,neumannCoeffs)
    enddo
    do i = 1,nEn
        call getNeumannIndices(iEn(i,:),ind,dEn(i),neumannLength)
        call applyNeumann(E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),dEn(i),neumannLength,neumannCoeffs)
    enddo
    do i = 1,nPn
        call getNeumannIndices(iPn(i,:),ind,dPn(i),neumannLength)
        
        P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = R(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))*E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
        call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),dPn(i),neumannLength,neumannCoeffs)
        R(iPn(i,1):iPn(i,2),iPn(i,3):iPn(i,4),iPn(i,5):iPn(i,6)) = P(iPn(i,1):iPn(i,2),iPn(i,3):iPn(i,4),iPn(i,5):iPn(i,6))/E(iPn(i,1):iPn(i,2),iPn(i,3):iPn(i,4),iPn(i,5):iPn(i,6))
    enddo
    
    ! CORNERS
    do i = 1,cN
        ! Temperature
        if(cAdiabatic(i).eq.1) then
            E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = 0
            select case(cD(i,1))
                case(1)
                    call getNeumannIndices(cL(i,:),ind,1,neumannLength)
                    P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) ! P is used just as a temporary variable here
                    call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),1,neumannLength,neumannCoeffs)
                    E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) + P(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))
                case(-1)
                    call getNeumannIndices(cL(i,:),ind,2,neumannLength)
                    P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
                    call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),2,neumannLength,neumannCoeffs)
                    E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) + P(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))
            end select
            
            select case(cD(i,2))
                case(1)
                    call getNeumannIndices(cL(i,:),ind,3,neumannLength)
                    P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
                    call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),3,neumannLength,neumannCoeffs)
                    E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) + P(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))
                case(-1)
                    call getNeumannIndices(cL(i,:),ind,4,neumannLength)
                    P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
                    call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),4,neumannLength,neumannCoeffs)
                    E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) + P(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))
            end select
            
            select case(cD(i,3))
                case(1)
                    call getNeumannIndices(cL(i,:),ind,5,neumannLength)
                    P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
                    call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),5,neumannLength,neumannCoeffs)
                    E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) + P(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))
                case(-1)
                    call getNeumannIndices(cL(i,:),ind,6,neumannLength)
                    P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
                    call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),6,neumannLength,neumannCoeffs)
                    E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) + P(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))
            end select
            E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))/sum(abs(cD(i,:)))
        endif
    
        ! Pressure
        R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = 0
        select case(cD(i,1))
            case(1)
                call getNeumannIndices(cL(i,:),ind,1,neumannLength)
                P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = R(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))*E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
                call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),1,neumannLength,neumannCoeffs)
                R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) + P(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))
            case(-1)
                call getNeumannIndices(cL(i,:),ind,2,neumannLength)
                P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = R(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))*E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
                call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),2,neumannLength,neumannCoeffs)
                R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) + P(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))
        end select
        
        select case(cD(i,2))
            case(1)
                call getNeumannIndices(cL(i,:),ind,3,neumannLength)
                P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = R(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))*E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
                call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),3,neumannLength,neumannCoeffs)
                R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) + P(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))
            case(-1)
                call getNeumannIndices(cL(i,:),ind,4,neumannLength)
                P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = R(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))*E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
                call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),4,neumannLength,neumannCoeffs)
                R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) + P(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))
        end select
        
        select case(cD(i,3))
            case(1)
                call getNeumannIndices(cL(i,:),ind,5,neumannLength)
                P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = R(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))*E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
                call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),5,neumannLength,neumannCoeffs)
                R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) + P(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))
            case(-1)
                call getNeumannIndices(cL(i,:),ind,6,neumannLength)
                P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = R(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))*E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
                call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),6,neumannLength,neumannCoeffs)
                R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) + P(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))
        end select
        R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6)) = R(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))/E(cL(i,1):cL(i,2),cL(i,3):cL(i,4),cL(i,5):cL(i,6))/sum(abs(cD(i,:)))
    enddo
    
    ! NULL SECOND DERIVATIVE CONDITIONS
    do i = 1,nUs
        call getNeumannIndices(iUs(i,:),ind,dUs(i),neumann2Length)
        call applyNeumann(U(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),dUs(i),neumann2Length,neumann2Coeffs)
    enddo
    do i = 1,nVs
        call getNeumannIndices(iVs(i,:),ind,dVs(i),neumann2Length)
        call applyNeumann(V(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),dVs(i),neumann2Length,neumann2Coeffs)
    enddo
    do i = 1,nWs
        call getNeumannIndices(iWs(i,:),ind,dWs(i),neumann2Length)
        call applyNeumann(W(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),dWs(i),neumann2Length,neumann2Coeffs)
    enddo
    do i = 1,nEs
        call getNeumannIndices(iEs(i,:),ind,dEs(i),neumann2Length)
        call applyNeumann(E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),dEs(i),neumann2Length,neumann2Coeffs)
    enddo
    do i = 1,nPs
        call getNeumannIndices(iPs(i,:),ind,dPs(i),neumann2Length)
        
        P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = R(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))*E(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
        call applyNeumann(P(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)),dPs(i),neumann2Length,neumann2Coeffs)
        R(iPs(i,1):iPs(i,2),iPs(i,3):iPs(i,4),iPs(i,5):iPs(i,6)) = P(iPs(i,1):iPs(i,2),iPs(i,3):iPs(i,4),iPs(i,5):iPs(i,6))/E(iPs(i,1):iPs(i,2),iPs(i,3):iPs(i,4),iPs(i,5):iPs(i,6))
    enddo
    
    ! CALL DISTURBANCES
    include 'runDisturbances.F90'
    
    end subroutine
    
    subroutine getNeumannIndices(indIn,ind,dir,neumannLength)

    implicit none
    
    integer, dimension(:), intent(in) :: indIn
    integer, dimension(6), intent(out) :: ind
    integer, intent(in) :: dir, neumannLength
    
    ind = indIn
    
    select case(dir)
        case (1)
            ind(2) = ind(2) + neumannLength
            
        case (2)
            ind(1) = ind(1) - neumannLength
        
        case (3)
            ind(4) = ind(4) + neumannLength
            
        case (4)
            ind(3) = ind(3) - neumannLength
            
        case (5)
            ind(6) = ind(6) + neumannLength
            
        case (6)
            ind(5) = ind(5) - neumannLength
    end select
    
    endsubroutine
    
    subroutine applyNeumann(U,dir,neumannLength,neumannCoeffs)
    
    implicit none
    
    real*8,  dimension(:,:,:), intent(inout) :: U
    integer, intent(in) :: dir, neumannLength
    real*8, dimension(:), intent(in) :: neumannCoeffs
    integer :: i
    
    select case(dir)
        case (1)
            U(1,:,:) = neumannCoeffs(1) * U(2,:,:)
            do i = 2,neumannLength
                U(1,:,:) = U(1,:,:) + neumannCoeffs(i) * U(i+1,:,:)
            enddo
            
        case (2)
            U(neumannLength+1,:,:) = neumannCoeffs(1) * U(neumannLength,:,:)
            do i = 2,neumannLength
                U(neumannLength+1,:,:) = U(neumannLength+1,:,:) + neumannCoeffs(i) * U(neumannLength-i+1,:,:)
            enddo
        
        case (3)
            U(:,1,:) = neumannCoeffs(1) * U(:,2,:)
            do i = 2,neumannLength
                U(:,1,:) = U(:,1,:) + neumannCoeffs(i) * U(:,i+1,:)
            enddo
            
        case (4)
            U(:,neumannLength+1,:) = neumannCoeffs(1) * U(:,neumannLength,:)
            do i = 2,neumannLength
                U(:,neumannLength+1,:) = U(:,neumannLength+1,:) + neumannCoeffs(i) * U(:,neumannLength-i+1,:)
            enddo
        
        case (5)
            U(:,:,1) = neumannCoeffs(1) * U(:,:,2)
            do i = 2,neumannLength
                U(:,:,1) = U(:,:,1) + neumannCoeffs(i) * U(:,:,i+1)
            enddo
            
        case (6)
            U(:,:,neumannLength+1) = neumannCoeffs(1) * U(:,:,neumannLength)
            do i = 2,neumannLength
                U(:,:,neumannLength+1) = U(:,:,neumannLength+1) + neumannCoeffs(i) * U(:,:,neumannLength-i+1)
            enddo
    end select
    
    endsubroutine
	
    end module
