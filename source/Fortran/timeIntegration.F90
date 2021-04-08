    module timeIntegration
    
    use decomp_2d
    use nseq
    
    contains

    subroutine RK4_2D(U,V,W,R,E,Unew,Vnew,Rnew,Enew, &
                & t, dt, Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                & derivnRHSx,derivnRHSy, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsDX, derivsDY, &
                & periodicX, periodicY, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    implicit none

    ! INPUTS
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: U,V,R,E
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(inout) :: W
    real*8,  intent(in)                   :: t, dt, Re, Ma, Pr, gamma, T0
    integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksY
    integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksY
    integer, intent(in)                   :: derivnRHSx, derivnRHSy
    real*8,  dimension(:,:), intent(in)   :: derivsAX, derivsBX, derivsCX, derivsAY, derivsBY, derivsCY
    real*8,  dimension(:,:), intent(in)   :: derivsDX, derivsDY
    real*8,  dimension(:,:,:), intent(in) :: derivsRX, derivsRY
    integer, intent(in)                   :: periodicX, periodicY
    integer, intent(in)                   :: neumannLength, neumann2Length
    real*8,  dimension(:), intent(in)     :: neumannCoeffs, neumann2Coeffs
    integer, intent(in)                   :: nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs
    integer, dimension(:,:), intent(in)   :: iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs
    real*8,  dimension(:), intent(in)     :: vUd, vVd, vWd, vPd, vEd
    integer, dimension(:), intent(in)     :: dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs
    integer, intent(in)                   :: cN
    integer, dimension(:,:), intent(in)   :: cL, cD
    integer, dimension(:), intent(in)     :: cAdiabatic
    
    ! OUTPUTS
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Unew,Vnew,Rnew,Enew
    
    ! INTERNAL VARIABLES
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: dU, dV, dR, dE
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: deltaU, deltaV, deltaR, deltaE
    
    ! Step 1
    call nseq2d(U,V,R,E,dU,dV,dR,dE, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                & derivnRHSx,derivnRHSy, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsDX, derivsDY, &
                & periodicX, periodicY)

    deltaU = dt/6.0d0*dU
    deltaV = dt/6.0d0*dV
    deltaR = dt/6.0d0*dR
    deltaE = dt/6.0d0*dE

    ! Step 2
    Unew = U + dt/2.0d0*dU
    Vnew = V + dt/2.0d0*dV
    Rnew = R + dt/2.0d0*dR
    Enew = E + dt/2.0d0*dE
    
    call applyBoundCond(Unew,Vnew,W,Rnew,Enew,gamma,t + dt/2.0, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)
                
    call nseq2d(Unew,Vnew,Rnew,Enew,dU,dV,dR,dE, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                & derivnRHSx,derivnRHSy, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsDX, derivsDY, &
                & periodicX, periodicY)

    deltaU = deltaU + dt/3.0d0*dU
    deltaV = deltaV + dt/3.0d0*dV
    deltaR = deltaR + dt/3.0d0*dR
    deltaE = deltaE + dt/3.0d0*dE
    
    ! Step 3
    Unew = U + dt/2.0d0*dU
    Vnew = V + dt/2.0d0*dV
    Rnew = R + dt/2.0d0*dR
    Enew = E + dt/2.0d0*dE

    call applyBoundCond(Unew,Vnew,W,Rnew,Enew,gamma, t + dt/2.0, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    call nseq2d(Unew,Vnew,Rnew,Enew,dU,dV,dR,dE, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                & derivnRHSx,derivnRHSy, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsDX, derivsDY, &
                & periodicX, periodicY)

    deltaU = deltaU + dt/3.0d0*dU
    deltaV = deltaV + dt/3.0d0*dV
    deltaR = deltaR + dt/3.0d0*dR
    deltaE = deltaE + dt/3.0d0*dE
    
    ! Step 4
    Unew = U + dt*dU
    Vnew = V + dt*dV
    Rnew = R + dt*dR
    Enew = E + dt*dE

    call applyBoundCond(Unew,Vnew,W,Rnew,Enew,gamma,t + dt, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    call nseq2d(Unew,Vnew,Rnew,Enew,dU,dV,dR,dE, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                & derivnRHSx,derivnRHSy, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsDX, derivsDY, &
                & periodicX, periodicY)
                
    Unew = U + deltaU + dt/6.0d0*dU
    Vnew = V + deltaV + dt/6.0d0*dV
    Rnew = R + deltaR + dt/6.0d0*dR
    Enew = E + deltaE + dt/6.0d0*dE
    
    ! Apply boundaries
    call applyBoundCond(Unew,Vnew,W,Rnew,Enew,gamma,t + dt, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    end subroutine
    
    subroutine EULER_2D(U,V,W,R,E,Unew,Vnew,Rnew,Enew, &
                & t, dt, Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                & derivnRHSx,derivnRHSy, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsDX, derivsDY, &
                & periodicX, periodicY, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    implicit none

    ! INPUTS
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: U,V,R,E
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(inout) :: W
    real*8,  intent(in)                   :: t, dt, Re, Ma, Pr, gamma, T0
    integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksY
    integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksY
    integer, intent(in)                   :: derivnRHSx, derivnRHSy
    real*8,  dimension(:,:), intent(in)   :: derivsAX, derivsBX, derivsCX, derivsAY, derivsBY, derivsCY
    real*8,  dimension(:,:), intent(in)   :: derivsDX, derivsDY
    real*8,  dimension(:,:,:), intent(in) :: derivsRX, derivsRY
    integer, intent(in)                   :: periodicX, periodicY
    integer, intent(in)                   :: neumannLength, neumann2Length
    real*8,  dimension(:), intent(in)     :: neumannCoeffs, neumann2Coeffs
    integer, intent(in)                   :: nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs
    integer, dimension(:,:), intent(in)   :: iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs
    real*8,  dimension(:), intent(in)     :: vUd, vVd, vWd, vPd, vEd
    integer, dimension(:), intent(in)     :: dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs
    integer, intent(in)                   :: cN
    integer, dimension(:,:), intent(in)   :: cL, cD
    integer, dimension(:), intent(in)     :: cAdiabatic
    
    ! OUTPUTS
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Unew,Vnew,Rnew,Enew
    
    ! INTERNAL VARIABLES
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: dU, dV, dR, dE
    
    call nseq2d(U,V,R,E,dU,dV,dR,dE, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                & derivnRHSx,derivnRHSy, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsDX, derivsDY, &
                & periodicX, periodicY)
                
    Unew = U + dt*dU
    Vnew = V + dt*dV
    Rnew = R + dt*dR
    Enew = E + dt*dE
    
    call applyBoundCond(Unew,Vnew,W,Rnew,Enew,gamma,t + dt, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)
    
    end subroutine EULER_2D
	
	
    subroutine SSPRK3_2D(U,V,W,R,E,Unew,Vnew,Rnew,Enew, &
                & t, dt, Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                & derivnRHSx,derivnRHSy, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsDX, derivsDY, &
                & periodicX, periodicY, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    implicit none

    ! INPUTS
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: U,V,R,E
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(inout) :: W
    real*8,  intent(in)                   :: t, dt, Re, Ma, Pr, gamma, T0
    integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksY
    integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksY
    integer, intent(in)                   :: derivnRHSx, derivnRHSy
    real*8,  dimension(:,:), intent(in)   :: derivsAX, derivsBX, derivsCX, derivsAY, derivsBY, derivsCY
    real*8,  dimension(:,:), intent(in)   :: derivsDX, derivsDY
    real*8,  dimension(:,:,:), intent(in) :: derivsRX, derivsRY
    integer, intent(in)                   :: periodicX, periodicY
    integer, intent(in)                   :: neumannLength, neumann2Length
    real*8,  dimension(:), intent(in)     :: neumannCoeffs, neumann2Coeffs
    integer, intent(in)                   :: nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs
    integer, dimension(:,:), intent(in)   :: iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs
    real*8,  dimension(:), intent(in)     :: vUd, vVd, vWd, vPd, vEd
    integer, dimension(:), intent(in)     :: dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs
    integer, intent(in)                   :: cN
    integer, dimension(:,:), intent(in)   :: cL, cD
    integer, dimension(:), intent(in)     :: cAdiabatic
    
    ! OUTPUTS
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Unew,Vnew,Rnew,Enew
    
    ! INTERNAL VARIABLES
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: Uk1, Vk1, Rk1, Ek1, Uk2, Vk2, Rk2, Ek2, Uk3, Vk3, Rk3, Ek3
    
	! Step 1
    call nseq2d(U,V,R,E,Uk1,Vk1,Rk1,Ek1, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                & derivnRHSx,derivnRHSy, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsDX, derivsDY, &
                & periodicX, periodicY)
                
    Unew = U + dt*Uk1
    Vnew = V + dt*Vk1
    Rnew = R + dt*Rk1
    Enew = E + dt*Ek1
    
    call applyBoundCond(Unew,Vnew,W,Rnew,Enew,gamma,t + dt, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)
				
	! Step 2
    call nseq2d(Unew,Vnew,Rnew,Enew,Uk2,Vk2,Rk2,Ek2, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                & derivnRHSx,derivnRHSy, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsDX, derivsDY, &
                & periodicX, periodicY)
                
    Unew = U + dt*(Uk1 + Uk2)/4.d0
    Vnew = V + dt*(Vk1 + Vk2)/4.d0
    Rnew = R + dt*(Rk1 + Rk2)/4.d0
    Enew = E + dt*(Ek1 + Ek2)/4.d0
    
    call applyBoundCond(Unew,Vnew,W,Rnew,Enew,gamma,t + dt/2.d0, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)
				
				
	! Step 3
    call nseq2d(Unew,Vnew,Rnew,Enew,Uk3,Vk3,Rk3,Ek3, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,derivBlocksX,derivBlocksY, &
                & derivnRHSx,derivnRHSy, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsDX, derivsDY, &
                & periodicX, periodicY)
                
    Unew = U + dt*(Uk1 + Uk2 + 4.d0*Uk3)/6.d0
    Vnew = V + dt*(Vk1 + Vk2 + 4.d0*Vk3)/6.d0
    Rnew = R + dt*(Rk1 + Rk2 + 4.d0*Rk3)/6.d0
    Enew = E + dt*(Ek1 + Ek2 + 4.d0*Ek3)/6.d0
    
    call applyBoundCond(Unew,Vnew,W,Rnew,Enew,gamma,t + dt, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)
    
    end subroutine SSPRK3_2D	
	
	
    subroutine RK4_3D(U,V,W,R,E,Unew,Vnew,Wnew,Rnew,Enew, &
                & t, dt, Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                & derivnRHSx,derivnRHSy,derivnRHSz, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                & derivsDX, derivsDY, derivsDZ, &
                & periodicX, periodicY, periodicZ, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    implicit none

    ! INPUTS
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: U,V,W,R,E
    real*8,  intent(in)                   :: t, dt, Re, Ma, Pr, gamma, T0
    integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ
    integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksY,derivBlocksZ
    integer, intent(in)                   :: derivnRHSx, derivnRHSy, derivnRHSz
    real*8,  dimension(:,:), intent(in)   :: derivsAX, derivsBX, derivsCX, derivsAY, derivsBY, derivsCY, derivsAZ, derivsBZ, derivsCZ
    real*8,  dimension(:,:), intent(in)   :: derivsDX, derivsDY, derivsDZ
    real*8,  dimension(:,:,:), intent(in) :: derivsRX, derivsRY, derivsRZ
    integer, intent(in)                   :: periodicX, periodicY, periodicZ
    integer, intent(in)                   :: neumannLength, neumann2Length
    real*8,  dimension(:), intent(in)     :: neumannCoeffs, neumann2Coeffs
    integer, intent(in)                   :: nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs
    integer, dimension(:,:), intent(in)   :: iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs
    real*8,  dimension(:), intent(in)     :: vUd, vVd, vWd, vPd, vEd
    integer, dimension(:), intent(in)     :: dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs
    integer, intent(in)                   :: cN
    integer, dimension(:,:), intent(in)   :: cL, cD
    integer, dimension(:), intent(in)     :: cAdiabatic
    
    ! OUTPUTS
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Unew,Vnew,Wnew,Rnew,Enew
    
    ! INTERNAL VARIABLES
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: dU, dV, dW, dR, dE
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: deltaU, deltaV, deltaW, deltaR, deltaE
    
    ! Step 1
    call nseq3d(U,V,W,R,E,dU,dV,dW,dR,dE, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                & derivnRHSx,derivnRHSy,derivnRHSz, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                & derivsDX, derivsDY, derivsDZ, &
                & periodicX, periodicY, periodicZ)

    deltaU = dt/6.0d0*dU
    deltaV = dt/6.0d0*dV
    deltaW = dt/6.0d0*dW
    deltaR = dt/6.0d0*dR
    deltaE = dt/6.0d0*dE

    ! Step 2
    Unew = U + dt/2.0d0*dU
    Vnew = V + dt/2.0d0*dV
    Wnew = W + dt/2.0d0*dW
    Rnew = R + dt/2.0d0*dR
    Enew = E + dt/2.0d0*dE
    
    call applyBoundCond(Unew,Vnew,Wnew,Rnew,Enew,gamma,t + dt/2.0, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)
                
    call nseq3d(Unew,Vnew,Wnew,Rnew,Enew,dU,dV,dW,dR,dE, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                & derivnRHSx,derivnRHSy,derivnRHSz, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                & derivsDX, derivsDY, derivsDZ, &
                & periodicX, periodicY, periodicZ)

    deltaU = deltaU + dt/3.0d0*dU
    deltaV = deltaV + dt/3.0d0*dV
    deltaW = deltaW + dt/3.0d0*dW
    deltaR = deltaR + dt/3.0d0*dR
    deltaE = deltaE + dt/3.0d0*dE
    
    ! Step 3
    Unew = U + dt/2.0d0*dU
    Vnew = V + dt/2.0d0*dV
    Wnew = W + dt/2.0d0*dW
    Rnew = R + dt/2.0d0*dR
    Enew = E + dt/2.0d0*dE

    call applyBoundCond(Unew,Vnew,Wnew,Rnew,Enew,gamma, t + dt/2.0, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    call nseq3d(Unew,Vnew,Wnew,Rnew,Enew,dU,dV,dW,dR,dE, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                & derivnRHSx,derivnRHSy,derivnRHSz, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                & derivsDX, derivsDY, derivsDZ, &
                & periodicX, periodicY, periodicZ)

    deltaU = deltaU + dt/3.0d0*dU
    deltaV = deltaV + dt/3.0d0*dV
    deltaW = deltaW + dt/3.0d0*dW
    deltaR = deltaR + dt/3.0d0*dR
    deltaE = deltaE + dt/3.0d0*dE
    
    ! Step 4
    Unew = U + dt*dU
    Vnew = V + dt*dV
    Wnew = W + dt*dW
    Rnew = R + dt*dR
    Enew = E + dt*dE

    call applyBoundCond(Unew,Vnew,Wnew,Rnew,Enew,gamma,t + dt, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    call nseq3d(Unew,Vnew,Wnew,Rnew,Enew,dU,dV,dW,dR,dE, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                & derivnRHSx,derivnRHSy,derivnRHSz, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                & derivsDX, derivsDY, derivsDZ, &
                & periodicX, periodicY, periodicZ)
                
    Unew = U + deltaU + dt/6.0d0*dU
    Vnew = V + deltaV + dt/6.0d0*dV
    Wnew = W + deltaW + dt/6.0d0*dW
    Rnew = R + deltaR + dt/6.0d0*dR
    Enew = E + deltaE + dt/6.0d0*dE

    
    ! Apply boundaries
    call applyBoundCond(Unew,Vnew,Wnew,Rnew,Enew,gamma,t + dt, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    end subroutine
    
    subroutine EULER_3D(U,V,W,R,E,Unew,Vnew,Wnew,Rnew,Enew, &
                & t, dt, Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                & derivnRHSx,derivnRHSy,derivnRHSz, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                & derivsDX, derivsDY, derivsDZ, &
                & periodicX, periodicY, periodicZ, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    implicit none

    ! INPUTS
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: U,V,W,R,E
    real*8,  intent(in)                   :: t, dt, Re, Ma, Pr, gamma, T0
    integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ
    integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksY,derivBlocksZ
    integer, intent(in)                   :: derivnRHSx, derivnRHSy, derivnRHSz
    real*8,  dimension(:,:), intent(in)   :: derivsAX, derivsBX, derivsCX, derivsAY, derivsBY, derivsCY, derivsAZ, derivsBZ, derivsCZ
    real*8,  dimension(:,:), intent(in)   :: derivsDX, derivsDY, derivsDZ
    real*8,  dimension(:,:,:), intent(in) :: derivsRX, derivsRY, derivsRZ
    integer, intent(in)                   :: periodicX, periodicY, periodicZ
    integer, intent(in)                   :: neumannLength, neumann2Length
    real*8,  dimension(:), intent(in)     :: neumannCoeffs, neumann2Coeffs
    integer, intent(in)                   :: nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs
    integer, dimension(:,:), intent(in)   :: iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs
    real*8,  dimension(:), intent(in)     :: vUd, vVd, vWd, vPd, vEd
    integer, dimension(:), intent(in)     :: dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs
    integer, intent(in)                   :: cN
    integer, dimension(:,:), intent(in)   :: cL, cD
    integer, dimension(:), intent(in)     :: cAdiabatic
    
    ! OUTPUTS
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Unew,Vnew,Wnew,Rnew,Enew
    
    ! INTERNAL VARIABLES
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: dU, dV, dW, dR, dE
    
    call nseq3d(U,V,W,R,E,dU,dV,dW,dR,dE, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                & derivnRHSx,derivnRHSy,derivnRHSz, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                & derivsDX, derivsDY, derivsDZ, &
                & periodicX, periodicY, periodicZ)

    Unew = U + dt*dU
    Vnew = V + dt*dV
    Wnew = W + dt*dW
    Rnew = R + dt*dR
    Enew = E + dt*dE
    
    call applyBoundCond(Unew,Vnew,Wnew,Rnew,Enew,gamma,t, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)
                
    end subroutine
	
	
    subroutine SSPRK3_3D(U,V,W,R,E,Unew,Vnew,Wnew,Rnew,Enew, &
                & t, dt, Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                & derivnRHSx,derivnRHSy,derivnRHSz, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                & derivsDX, derivsDY, derivsDZ, &
                & periodicX, periodicY, periodicZ, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)

    implicit none

    ! INPUTS
    real*8,  dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: U,V,W,R,E
    real*8,  intent(in)                   :: t, dt, Re, Ma, Pr, gamma, T0
    integer, intent(in)                   :: nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ
    integer, dimension(:,:), intent(in)   :: derivBlocksX,derivBlocksY,derivBlocksZ
    integer, intent(in)                   :: derivnRHSx, derivnRHSy, derivnRHSz
    real*8,  dimension(:,:), intent(in)   :: derivsAX, derivsBX, derivsCX, derivsAY, derivsBY, derivsCY, derivsAZ, derivsBZ, derivsCZ
    real*8,  dimension(:,:), intent(in)   :: derivsDX, derivsDY, derivsDZ
    real*8,  dimension(:,:,:), intent(in) :: derivsRX, derivsRY, derivsRZ
    integer, intent(in)                   :: periodicX, periodicY, periodicZ
    integer, intent(in)                   :: neumannLength, neumann2Length
    real*8,  dimension(:), intent(in)     :: neumannCoeffs, neumann2Coeffs
    integer, intent(in)                   :: nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs
    integer, dimension(:,:), intent(in)   :: iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs
    real*8,  dimension(:), intent(in)     :: vUd, vVd, vWd, vPd, vEd
    integer, dimension(:), intent(in)     :: dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs
    integer, intent(in)                   :: cN
    integer, dimension(:,:), intent(in)   :: cL, cD
    integer, dimension(:), intent(in)     :: cAdiabatic
    
    ! OUTPUTS
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(out) :: Unew,Vnew,Wnew,Rnew,Enew
    
    ! INTERNAL VARIABLES
    real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: Uk1, Vk1, Wk1, Rk1, Ek1, Uk2, Vk2, Wk2, Rk2, Ek2, Uk3, Vk3, Wk3, Rk3, Ek3
    
	! Step 1
    call nseq3d(U,V,W,R,E,Uk1,Vk1,Wk1,Rk1,Ek1, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                & derivnRHSx,derivnRHSy,derivnRHSz, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                & derivsDX, derivsDY, derivsDZ, &
                & periodicX, periodicY, periodicZ)
                
    Unew = U + dt*Uk1
    Vnew = V + dt*Vk1
	Wnew = W + dt*Wk1
    Rnew = R + dt*Rk1
    Enew = E + dt*Ek1
    
    call applyBoundCond(Unew,Vnew,Wnew,Rnew,Enew,gamma,t + dt, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)
				
	! Step 2
    call nseq3d(Unew,Vnew,Wnew,Rnew,Enew,Uk2,Vk2,Wk2,Rk2,Ek2, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                & derivnRHSx,derivnRHSy,derivnRHSz, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                & derivsDX, derivsDY, derivsDZ, &
                & periodicX, periodicY, periodicZ)
                
    Unew = U + dt*(Uk1 + Uk2)/4.d0
    Vnew = V + dt*(Vk1 + Vk2)/4.d0
	Wnew = W + dt*(Wk1 + Wk2)/4.d0
    Rnew = R + dt*(Rk1 + Rk2)/4.d0
    Enew = E + dt*(Ek1 + Ek2)/4.d0
	
    call applyBoundCond(Unew,Vnew,Wnew,Rnew,Enew,gamma,t + dt/2.d0, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)
				
				
	! Step 3
    call nseq3d(Unew,Vnew,Wnew,Rnew,Enew,Uk3,Vk3,Wk3,Rk3,Ek3, &
                & Re, Ma, Pr, gamma, T0, &
                & nDerivBlocksX,nDerivBlocksY,nDerivBlocksZ,derivBlocksX,derivBlocksY,derivBlocksZ, &
                & derivnRHSx,derivnRHSy,derivnRHSz, &
                & derivsAX, derivsBX, derivsCX, derivsRX, &
                & derivsAY, derivsBY, derivsCY, derivsRY, &
                & derivsAZ, derivsBZ, derivsCZ, derivsRZ, &
                & derivsDX, derivsDY, derivsDZ, &
                & periodicX, periodicY, periodicZ)
                
    Unew = U + dt*(Uk1 + Uk2 + 4.d0*Uk3)/6.d0
    Vnew = V + dt*(Vk1 + Vk2 + 4.d0*Vk3)/6.d0
	Wnew = W + dt*(Wk1 + Wk2 + 4.d0*Wk3)/6.d0
    Rnew = R + dt*(Rk1 + Rk2 + 4.d0*Rk3)/6.d0
    Enew = E + dt*(Ek1 + Ek2 + 4.d0*Ek3)/6.d0
    
    call applyBoundCond(Unew,Vnew,Wnew,Rnew,Enew,gamma,t + dt, &
                & neumannLength, neumann2Length, neumannCoeffs, neumann2Coeffs, &
                & nUd, nVd, nWd, nPd, nEd, nUn, nVn, nWn, nPn, nEn, nUs, nVs, nWs, nPs, nEs, &
                & iUd, iVd, iWd, iPd, iEd, iUn, iVn, iWn, iPn, iEn, iUs, iVs, iWs, iPs, iEs, &
                & vUd, vVd, vWd, vPd, vEd, dUn, dVn, dWn, dPn, dEn, dUs, dVs, dWs, dPs, dEs, &
                & cN, cL, cD, cAdiabatic)
    
    end subroutine SSPRK3_3D
	
	end module