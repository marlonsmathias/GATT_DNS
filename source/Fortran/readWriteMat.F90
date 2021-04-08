#include "fintrf.h"

module readWriteMat

  contains

    subroutine readFlow(timeStep,t,U,V,W,R,E)

    ! DECLARE VARIABLES
    implicit none

    ! Inputs
    integer, intent(in) :: timeStep

    ! Outputs
    real*8, intent(out) :: t
    real*8, dimension(:,:,:), allocatable, intent(out) :: U, V, W, R, E

    ! Internals
    mwPointer :: matOpen, matClose, mxGetData, matGetVariable
    mwPointer :: filePointer, infoPointer
    mwPointer :: varPointerT, dataPointerT
    mwPointer :: varPointerU, dataPointerU
    mwPointer :: varPointerV, dataPointerV
    mwPointer :: varPointerW, dataPointerW
    mwPointer :: varPointerR, dataPointerR
    mwPointer :: varPointerE, dataPointerE
    mwPointer :: matGetVariableInfo, mxGetDimensions
    mwPointer :: mxGetNumberOfDimensions
    mwSize :: varSize, varSize1
    mwSize :: numberOfDimensions
    character(len=10) :: timeChar

    integer*8 :: dims(3)
    
    ! OPEN MAT FILE
    write(timeChar,fmt='(i10.10)') timeStep
    filePointer = matOpen('../flow_'//timeChar//'.mat', 'r')

    ! GET POINTERS
    varPointerT = matGetVariable(filePointer, 't')
    dataPointerT = mxGetData(varPointerT)
    varPointerU = matGetVariable(filePointer, 'U')
    dataPointerU = mxGetData(varPointerU)
    varPointerV = matGetVariable(filePointer, 'V')
    dataPointerV = mxGetData(varPointerV)
    varPointerW = matGetVariable(filePointer, 'W')
    dataPointerW = mxGetData(varPointerW)
    varPointerR = matGetVariable(filePointer, 'R')
    dataPointerR = mxGetData(varPointerR)
    varPointerE = matGetVariable(filePointer, 'E')
    dataPointerE = mxGetData(varPointerE)
    
    ! GET DIMENSIONS
    infoPointer = mxGetDimensions(varPointerU)
    numberOfDimensions = mxGetNumberOfDimensions(varPointerU)
    call mxCopyPtrToInteger8(infoPointer, dims, numberOfDimensions)

    ! If file has only two dimensions, set the third to one
    if(numberOfDimensions.eq.2) then
        dims(3) = 1
    end if
    
    varSize = dims(1)*dims(2)*dims(3)
    varSize1 = 1

    ! ALLOCATE OUTPUT VARIABLE
    allocate(U(dims(1),dims(2),dims(3)))
    allocate(V(dims(1),dims(2),dims(3)))
    allocate(W(dims(1),dims(2),dims(3)))
    allocate(R(dims(1),dims(2),dims(3)))
    allocate(E(dims(1),dims(2),dims(3)))

    ! READ VALUES
    call mxCopyPtrToReal8(dataPointerT,t,varSize1)
    call mxCopyPtrToReal8(dataPointerU,U,varSize)
    call mxCopyPtrToReal8(dataPointerV,V,varSize)
    call mxCopyPtrToReal8(dataPointerW,W,varSize)
    call mxCopyPtrToReal8(dataPointerR,R,varSize)
    call mxCopyPtrToReal8(dataPointerE,E,varSize)

    ! CLOSE MAT FILE
    filePointer = matClose(filePointer)
    
    end subroutine
    
    subroutine readSFD(sfd_X)
    
    ! DECLARE VARIABLES
    implicit none

    ! Outputs
    real*8, dimension(:,:,:), allocatable, intent(out) :: sfd_X

    ! Internals
    mwPointer :: matOpen, matClose, mxGetData, matGetVariable
    mwPointer :: filePointer, infoPointer
    mwPointer :: varPointer, dataPointer
    mwPointer :: matGetVariableInfo, mxGetDimensions
    mwPointer :: mxGetNumberOfDimensions
    mwSize :: varSize
    mwSize :: numberOfDimensions

    integer*8 :: dims(3)
    
    ! OPEN MAT FILE
    
    filePointer = matOpen('SFD.mat', 'r')

    ! GET POINTERS
    varPointer = matGetVariable(filePointer, 'SFD_X')
    dataPointer = mxGetData(varPointer)
    
    ! GET DIMENSIONS
    infoPointer = mxGetDimensions(varPointer)
    numberOfDimensions = mxGetNumberOfDimensions(varPointer)
    call mxCopyPtrToInteger8(infoPointer, dims, numberOfDimensions)

    ! If file has only two dimensions, set the third to one
    if(numberOfDimensions.eq.2) then
        dims(3) = 1
    end if
    
    varSize = dims(1)*dims(2)*dims(3)

    ! ALLOCATE OUTPUT VARIABLE
    allocate(sfd_X(dims(1),dims(2),dims(3)))

    ! READ VALUES
    call mxCopyPtrToReal8(dataPointer,sfd_X,varSize)

    ! CLOSE MAT FILE
    filePointer = matClose(filePointer)
    
    end subroutine
    
    
    subroutine writeFlow(timeStep,t,U,V,W,R,E)


    ! DECLARE VARIABLES
    implicit none

    ! Inputs
    integer, intent(in) :: timeStep
    real*8, intent(in) :: t
    real*8, dimension(:,:,:), intent(in) :: U, V, W, R, E

    ! Internals
    mwPointer :: matOpen, matClose, matPutVariable, mxGetPr
    mwPointer :: filePointer, mxCreateNumericArray
    mwPointer :: varPointerT, dataPointerT
    mwPointer :: varPointerU, dataPointerU
    mwPointer :: varPointerV, dataPointerV
    mwPointer :: varPointerW, dataPointerW
    mwPointer :: varPointerR, dataPointerR
    mwPointer :: varPointerE, dataPointerE
    mwSize :: fileSize, fileSize1
    mwSize :: numberOfDimensions, numberOfDimensions1
    integer*8 :: dims(3), dims1
    integer*4 :: writeStatus, classID, mxClassIDFromClassName, complexFlag
    character(len=10) :: timeChar

    ! GET DATA SIZE
    numberOfDimensions = 3
    numberOfDimensions1 = 1
    dims = shape(U)
    dims1 = 1
    fileSize = dims(1)*dims(2)*dims(3)
    fileSize1 = 1

    ! OPEN MAT FILE
    write(timeChar,fmt='(i10.10)') timeStep
	
	if(timeStep.lt.0) then
		filePointer = matOpen('../flow_diverged.mat', 'w7.3')
	else
		filePointer = matOpen('../flow_'//timeChar//'.mat', 'w7.3')
	end if
	
    

    ! GET POINTERS
    classID = mxClassIDFromClassName('double')
    complexFlag = 0
    
    dataPointerT = mxCreateNumericArray(numberOfDimensions1, dims1, classID, complexFlag)
    varPointerT = mxGetPr(dataPointerT)
    dataPointerU = mxCreateNumericArray(numberOfDimensions, dims, classID, complexFlag)
    varPointerU = mxGetPr(dataPointerU)
    dataPointerV = mxCreateNumericArray(numberOfDimensions, dims, classID, complexFlag)
    varPointerV = mxGetPr(dataPointerV)
    dataPointerW = mxCreateNumericArray(numberOfDimensions, dims, classID, complexFlag)
    varPointerW = mxGetPr(dataPointerW)
    dataPointerR = mxCreateNumericArray(numberOfDimensions, dims, classID, complexFlag)
    varPointerR = mxGetPr(dataPointerR)
    dataPointerE = mxCreateNumericArray(numberOfDimensions, dims, classID, complexFlag)
    varPointerE = mxGetPr(dataPointerE)

    ! COPY DATA TO POINTER
    call mxCopyReal8ToPtr(t, varPointerT, fileSize1)
    call mxCopyReal8ToPtr(U, varPointerU, fileSize)
    call mxCopyReal8ToPtr(V, varPointerV, fileSize)
    call mxCopyReal8ToPtr(W, varPointerW, fileSize)
    call mxCopyReal8ToPtr(R, varPointerR, fileSize)
    call mxCopyReal8ToPtr(E, varPointerE, fileSize)

    ! WRITE VALUES
    writeStatus = matPutVariable(filePointer, 't', dataPointerT)
    writeStatus = matPutVariable(filePointer, 'U', dataPointerU)
    writeStatus = matPutVariable(filePointer, 'V', dataPointerV)
    writeStatus = matPutVariable(filePointer, 'W', dataPointerW)
    writeStatus = matPutVariable(filePointer, 'R', dataPointerR)
    writeStatus = matPutVariable(filePointer, 'E', dataPointerE)

	! DEALLOCATE MEMORY
	call mxDestroyArray(dataPointerT)
	call mxDestroyArray(dataPointerU)
	call mxDestroyArray(dataPointerV)
	call mxDestroyArray(dataPointerW)
	call mxDestroyArray(dataPointerR)
	call mxDestroyArray(dataPointerE)

    ! CLOSE MAT FILE
    filePointer = matClose(filePointer)


    end subroutine
	
	
	subroutine readMeanFlow(U,V,W,R,E)

    ! DECLARE VARIABLES
    implicit none

    ! Outputs
    real*8, dimension(:,:,:), allocatable, intent(out) :: U, V, W, R, E

    ! Internals
    mwPointer :: matOpen, matClose, mxGetData, matGetVariable
    mwPointer :: filePointer, infoPointer
    mwPointer :: varPointerU, dataPointerU
    mwPointer :: varPointerV, dataPointerV
    mwPointer :: varPointerW, dataPointerW
    mwPointer :: varPointerR, dataPointerR
    mwPointer :: varPointerE, dataPointerE
    mwPointer :: matGetVariableInfo, mxGetDimensions
    mwPointer :: mxGetNumberOfDimensions
    mwSize :: varSize, varSize1
    mwSize :: numberOfDimensions

    integer*8 :: dims(3)
    
    ! OPEN MAT FILE
    filePointer = matOpen('../meanflowSFD.mat', 'r')

    ! GET POINTERS
    varPointerU = matGetVariable(filePointer, 'U')
    dataPointerU = mxGetData(varPointerU)
    varPointerV = matGetVariable(filePointer, 'V')
    dataPointerV = mxGetData(varPointerV)
    varPointerW = matGetVariable(filePointer, 'W')
    dataPointerW = mxGetData(varPointerW)
    varPointerR = matGetVariable(filePointer, 'R')
    dataPointerR = mxGetData(varPointerR)
    varPointerE = matGetVariable(filePointer, 'E')
    dataPointerE = mxGetData(varPointerE)
    
    ! GET DIMENSIONS
    infoPointer = mxGetDimensions(varPointerU)
    numberOfDimensions = mxGetNumberOfDimensions(varPointerU)
    call mxCopyPtrToInteger8(infoPointer, dims, numberOfDimensions)

    ! If file has only two dimensions, set the third to one
    if(numberOfDimensions.eq.2) then
        dims(3) = 1
    end if
    
    varSize = dims(1)*dims(2)*dims(3)
    varSize1 = 1

    ! ALLOCATE OUTPUT VARIABLE
    allocate(U(dims(1),dims(2),dims(3)))
    allocate(V(dims(1),dims(2),dims(3)))
    allocate(W(dims(1),dims(2),dims(3)))
    allocate(R(dims(1),dims(2),dims(3)))
    allocate(E(dims(1),dims(2),dims(3)))

    ! READ VALUES
    call mxCopyPtrToReal8(dataPointerU,U,varSize)
    call mxCopyPtrToReal8(dataPointerV,V,varSize)
    call mxCopyPtrToReal8(dataPointerW,W,varSize)
    call mxCopyPtrToReal8(dataPointerR,R,varSize)
    call mxCopyPtrToReal8(dataPointerE,E,varSize)

    ! CLOSE MAT FILE
    filePointer = matClose(filePointer)
    
    end subroutine
	
	
    subroutine writeMeanFlow(U,V,W,R,E)


    ! DECLARE VARIABLES
    implicit none

    ! Inputs
    real*8, dimension(:,:,:), intent(in) :: U, V, W, R, E

    ! Internals
    mwPointer :: matOpen, matClose, matPutVariable, mxGetPr
    mwPointer :: filePointer, mxCreateNumericArray
    mwPointer :: varPointerU, dataPointerU
    mwPointer :: varPointerV, dataPointerV
    mwPointer :: varPointerW, dataPointerW
    mwPointer :: varPointerR, dataPointerR
    mwPointer :: varPointerE, dataPointerE
    mwSize :: fileSize
    mwSize :: numberOfDimensions
    integer*8 :: dims(3)
    integer*4 :: writeStatus, classID, mxClassIDFromClassName, complexFlag

    ! GET DATA SIZE
    numberOfDimensions = 3
    dims = shape(U)
    fileSize = dims(1)*dims(2)*dims(3)

    ! OPEN MAT FILE
    filePointer = matOpen('../meanflowSFD.mat', 'w7.3')

    ! GET POINTERS
    classID = mxClassIDFromClassName('double')
    complexFlag = 0
    
    dataPointerU = mxCreateNumericArray(numberOfDimensions, dims, classID, complexFlag)
    varPointerU = mxGetPr(dataPointerU)
    dataPointerV = mxCreateNumericArray(numberOfDimensions, dims, classID, complexFlag)
    varPointerV = mxGetPr(dataPointerV)
    dataPointerW = mxCreateNumericArray(numberOfDimensions, dims, classID, complexFlag)
    varPointerW = mxGetPr(dataPointerW)
    dataPointerR = mxCreateNumericArray(numberOfDimensions, dims, classID, complexFlag)
    varPointerR = mxGetPr(dataPointerR)
    dataPointerE = mxCreateNumericArray(numberOfDimensions, dims, classID, complexFlag)
    varPointerE = mxGetPr(dataPointerE)

    ! COPY DATA TO POINTER
    call mxCopyReal8ToPtr(U, varPointerU, fileSize)
    call mxCopyReal8ToPtr(V, varPointerV, fileSize)
    call mxCopyReal8ToPtr(W, varPointerW, fileSize)
    call mxCopyReal8ToPtr(R, varPointerR, fileSize)
    call mxCopyReal8ToPtr(E, varPointerE, fileSize)

    ! WRITE VALUES
    writeStatus = matPutVariable(filePointer, 'U', dataPointerU)
    writeStatus = matPutVariable(filePointer, 'V', dataPointerV)
    writeStatus = matPutVariable(filePointer, 'W', dataPointerW)
    writeStatus = matPutVariable(filePointer, 'R', dataPointerR)
    writeStatus = matPutVariable(filePointer, 'E', dataPointerE)

	! DEALLOCATE MEMORY
	call mxDestroyArray(dataPointerU)
	call mxDestroyArray(dataPointerV)
	call mxDestroyArray(dataPointerW)
	call mxDestroyArray(dataPointerR)
	call mxDestroyArray(dataPointerE)

    ! CLOSE MAT FILE
    filePointer = matClose(filePointer)


    end subroutine

end module
