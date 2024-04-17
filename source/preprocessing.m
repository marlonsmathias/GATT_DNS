% This script is the main preprocessing routine, which will generate the mesh, prepare the boundaries, the matrices and write the Fortran files

%% Generate mesh
% Add fixed points to mesh structure
% The mesh will be slightly transformed so that key points are present at
% the correct positions if mesh.x.matchFixed = true
% X = 0 and Y = 0 are always added
meshAddFixedPoints

% Run mesh generator
[mesh.X, mesh.x, mesh.nx] = generateMesh(domain.xi,domain.xf,mesh.x,'X');
[mesh.Y, mesh.y, mesh.ny] = generateMesh(domain.yi,domain.yf,mesh.y,'Y');
[mesh.Z, mesh.z, mesh.nz] = generateMesh(domain.zi,domain.zf,mesh.z,'Z');

%% Select boundary conditions
[boundary,mesh] = getBoundaryConditions(flowType,mesh,flowParameters,[numMethods.neumannOrder numMethods.neumann2Order]);

domainSlicesY = getDomainSlices(mesh.ny,p_row);
domainSlicesZ = getDomainSlices(mesh.nz,p_col);

boundaryInfo = initBoundaries(boundary,mesh,domainSlicesY,domainSlicesZ,p_row,p_col);

%% Make matrices for spatial derivatives and for spatial filters
matrices = makeMatrices(mesh,domain,boundary,numMethods);

matrices.x = prepareThomas(matrices.x);
matrices.x.blocks = getMatrixTypeBlocks(matrices.x.types,p_row,p_col);
matrices.y = prepareThomas(matrices.y);
matrices.y.blocks = getMatrixTypeBlocks(matrices.y.types,p_row,p_col);
if mesh.nz>1
    matrices.z = prepareThomas(matrices.z);
    matrices.z.blocks = getMatrixTypeBlocks(matrices.z.types,p_row,p_col);
end

matrices.neumannCoeffs = boundary.neumannCoeffs;
matrices.neumann2Coeffs = boundary.neumann2Coeffs;

%% Write files

if ~exist(caseName,'dir')
    mkdir(caseName);
end

%% Prepare SFD
if isfield(numMethods,'SFD')
	if numMethods.SFD.type == 2
		calcSFDregion
	end
	if isinf(numMethods.SFD.Delta)
		numMethods.SFD.Delta = -1;
	end
	if numMethods.SFD.type > 0 && (exist([caseName '/meanflowSFD.mat'],'file') || isfield(flowType.initial,'meanFile'))
		if ~isfield(numMethods.SFD,'resume')
			numMethods.SFD.resume = 1;
		end
	else
		numMethods.SFD.resume = 0;
	end
end

% Mesh file
X = mesh.X; %#ok<*NASGU>
Y = mesh.Y;
Z = mesh.Z;
wall = boundary.insideWall;
save([caseName '/mesh.mat'], 'X', 'Y', 'Z', 'wall', 'flowParameters', 'flowType');
clear X Y Z

% Check for previous save files
if ~exist('runningLST','var')
	[nStep, nx, ny, nz] = checkPreviousRun(caseName);
	if ~isempty(nStep)
		genInitialFlow = false;
		time.nStep = nStep;
		if (nx ~= mesh.nx || ny ~= mesh.ny || nz ~= mesh.nz)
			error('Mesh size has changed since last run: Parameters file indicates %dx%dx%d but mat file contains data for %dx%dx%d', mesh.nx,mesh.ny,mesh.nz, nx,ny,nz)
		end
	else
		genInitialFlow = true;
		time.nStep = 0;
	end
	clear nx ny nz nStep
else
	genInitialFlow = false;
	time.nStep = 0;
	numMethods.SFD.resume = 0;
end

% Fortran files
if ~exist([caseName '/bin'],'dir')
    mkdir([caseName '/bin']);
end

disturbTypes = writeFortranDisturbances(caseName,boundaryInfo,tridimensional);

writeFortranParameters(caseName,mesh,flowParameters,time,numMethods,logAll,p_row,p_col)

writeFortranMatrices(caseName,matrices,numMethods,mesh)

writeFortranBoundaries(caseName,boundaryInfo)

% SFD
if numMethods.SFD.type == 2
    save([caseName '/bin/SFD.mat'], 'SFD_X');
end
