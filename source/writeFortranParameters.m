function writeFortranParameters(caseName,mesh,flowParameters,time,numMethods,logAll,p_row,p_col)

% This function writes the parameters.F90 file, which contains basic parameters that will be used in the simulation

    outFile = fopen([caseName '/bin/parameters.F90'],'w');
    
    fprintf(outFile,'    integer :: nx = %d\n',mesh.nx);
    fprintf(outFile,'    integer :: ny = %d\n',mesh.ny);
    fprintf(outFile,'    integer :: nz = %d\n\n',mesh.nz);
    
    fprintf(outFile,'    real*8 :: Re = %.20fd0\n',flowParameters.Re);
    fprintf(outFile,'    real*8 :: Ma = %.20fd0\n',flowParameters.Ma);
    fprintf(outFile,'    real*8 :: Pr = %.20fd0\n',flowParameters.Pr);
    fprintf(outFile,'    real*8 :: T0 = %.20fd0\n',flowParameters.T0);
    fprintf(outFile,'    real*8 :: gamma = %.20fd0\n\n',flowParameters.gamma);
    
    fprintf(outFile,'    real*8 :: dtmax = %.20fd0\n',time.dt);
    fprintf(outFile,'    real*8 :: maxCFL = %.20fd0\n',time.maxCFL);
    
	if logAll == 0
		logAll = 2147483647;
	end
    fprintf(outFile,'    integer :: logAll = %d\n', logAll);
	
	fprintf(outFile,'    integer :: nSave = %d\n', time.nStep);
    
	if ~isfield(mesh,'trackedNorm') || ~mesh.trackedNorm
		fprintf(outFile,'    real*8 :: trackedNorm = 0.d0\n');
	else
		fprintf(outFile,'    real*8 :: trackedNorm = %.20fd0\n',1/((flowParameters.gamma^2-flowParameters.gamma)*flowParameters.Ma^2));
	end
	
    switch time.control
        case 'dt'
            fprintf(outFile,'    integer :: timeControl = 1\n');
            fprintf(outFile,'    integer :: qTimesInt = %d\n',time.qtimes);
            fprintf(outFile,'    real*8  :: qTimesReal\n\n');
            fprintf(outFile,'    integer :: tmaxInt = %d\n',time.tmax);
            fprintf(outFile,'    real*8  :: tmaxReal\n\n');
        case 'cfl'
            fprintf(outFile,'    integer :: timeControl = 2\n');
            fprintf(outFile,'    integer :: qTimesInt\n');
            fprintf(outFile,'    real*8  :: qtimesReal = %.20fd0\n\n',time.qtimes);
            fprintf(outFile,'    integer :: tmaxInt\n');
            fprintf(outFile,'    real*8  :: tmaxReal = %.20fd0\n\n',time.tmax);
        otherwise
            error('Unrecognized type of time control. Use either dt or cfl')
    end
    
    switch numMethods.timeStepping
        case 'RK4'
            fprintf(outFile,'    integer :: timeStepping = 1\n');
        case 'Euler'
            fprintf(outFile,'    integer :: timeStepping = 2\n');
		case 'SSPRK3'
            fprintf(outFile,'    integer :: timeStepping = 3\n');
        otherwise
            error('Unrecognized time stepping method')
    end
    
	if isfield(numMethods,'SFD')
		fprintf(outFile,'    integer :: SFD = %d\n',numMethods.SFD.type);
		fprintf(outFile,'    real*8 :: SFD_Delta = %.20fd0\n',numMethods.SFD.Delta);
		fprintf(outFile,'    real*8 :: SFD_X_val = %.20fd0\n',numMethods.SFD.X);
		fprintf(outFile,'    integer :: resumeMeanFlow = %d\n\n',numMethods.SFD.resume);
	else
		fprintf(outFile,'    integer :: SFD = 0\n');
		fprintf(outFile,'    real*8 :: SFD_Delta = %.20fd0\n',0);
		fprintf(outFile,'    real*8 :: SFD_X_val = %.20fd0\n',0);
		fprintf(outFile,'    integer :: resumeMeanFlow = 0\n\n');
	end
	
	if ~isfield(numMethods,'spatialFilterTime') || numMethods.spatialFilterTime <= 0
		fprintf(outFile,'    real*8 :: FilterCharTime = %.20fd0\n',-1);
	else
		fprintf(outFile,'    real*8 :: FilterCharTime = %.20fd0\n',numMethods.spatialFilterTime);
	end
    
    if mesh.nz == 1
        dxmin = [1/min(diff(mesh.X)), 1/min(diff(mesh.Y)), 0];
    else
        dxmin = [1/min(diff(mesh.X)), 1/min(diff(mesh.Y)), 1/min(diff(mesh.Z))];
    end
	
	if isfield(time,'CFLignoreZ') && time.CFLignoreZ
		dxmin(3) = 0;
	end
	
    fprintf(outFile,'    real*8,dimension(3) :: dxmin = (/%.20fd0,%.20fd0,%.20fd0/)\n\n', dxmin(1), dxmin(2), dxmin(3));
    
    fprintf(outFile,'    integer :: p_row = %d\n',p_row);
    fprintf(outFile,'    integer :: p_col = %d\n',p_col);
    
    fclose(outFile);

end
