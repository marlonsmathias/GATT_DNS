function writeFortranMatrices(caseName,matrices,numMethods,mesh)

% This file writes matrices.F90, which contains all the matrices for derivatives and filters, as well as th regions to which they will be applied
% The coefficients for Neumann boundary conditions are also defined here

outFile = fopen([caseName '/bin/matrices.F90'],'w');

I = size(matrices.y.types,1);
J = size(matrices.x.types,1);
K = size(matrices.x.types,2);
nProcs = length(matrices.x.blocks);

%% Write blocks
fprintf(outFile,'    select case(nrank)\n');

for i = 1:nProcs
    fprintf(outFile,'        case(%d)\n',i-1);
    
    % For X
    blocks = matrices.x.blocks{i};
    nBlocks = size(blocks,1);
    fprintf(outFile,'            nDerivBlocksX = %d\n',nBlocks);
    fprintf(outFile,'            allocate(derivBlocksX(%d,5))\n',nBlocks);
    
    fprintf(outFile,'            derivBlocksX = reshape((/');
    for n = 1:nBlocks*5-1
        fprintf(outFile, '%d,', blocks(n));
    end
    fprintf(outFile, '%d', blocks(end));
    fprintf(outFile,'/),shape(derivBlocksX))\n\n');
    
    % For Y
    blocks = matrices.y.blocks{i};
    nBlocks = size(blocks,1);
    fprintf(outFile,'            nDerivBlocksY = %d\n',nBlocks);
    fprintf(outFile,'            allocate(derivBlocksY(%d,5))\n',nBlocks);
    
    fprintf(outFile,'        	 derivBlocksY = reshape((/');
    for n = 1:nBlocks*5-1
        fprintf(outFile, '%d,', blocks(n));
    end
    fprintf(outFile, '%d', blocks(end));
    fprintf(outFile,'/),shape(derivBlocksY))\n\n');
    
    % For Z
    if K>1
        blocks = matrices.z.blocks{i};
        nBlocks = size(blocks,1);
        fprintf(outFile,'            nDerivBlocksZ = %d\n',nBlocks);
        fprintf(outFile,'            allocate(derivBlocksZ(%d,5))\n',nBlocks);

        fprintf(outFile,'            derivBlocksZ = reshape((/');
        for n = 1:nBlocks*5-1
            fprintf(outFile, '%d,', blocks(n));
        end
        fprintf(outFile, '%d', blocks(end));
        fprintf(outFile,'/),shape(derivBlocksZ))\n\n');
    end
end

fprintf(outFile,'    end select\n\n');

%% Write filter info
fprintf(outFile,'    filterX = %d\n',numMethods.filterDirections(1));
fprintf(outFile,'    filterY = %d\n',numMethods.filterDirections(2));
fprintf(outFile,'    filterZ = %d\n',numMethods.filterDirections(3));

%% Write matrices
%% For X derivatives
fprintf(outFile,'    derivnRHSx = %d\n',matrices.x.nRHS);
fprintf(outFile,'    filternRHSx = %d\n',matrices.x.nRHSf);

fprintf(outFile,'    allocate(derivsAX(%d,%d))\n',I-1,matrices.x.nTypes);
fprintf(outFile,'    allocate(derivsBX(%d,%d))\n',I-matrices.x.periodic,matrices.x.nTypes);
fprintf(outFile,'    allocate(derivsCX(%d,%d))\n',I-1,matrices.x.nTypes);
fprintf(outFile,'    allocate(derivsRX(%d,%d,%d))\n\n',I,2*matrices.x.nRHS-1,matrices.x.nTypes);

fprintf(outFile,'    allocate(filterAX(%d,%d))\n',I-1,matrices.x.nTypes);
fprintf(outFile,'    allocate(filterBX(%d,%d))\n',I-matrices.x.periodic,matrices.x.nTypes);
fprintf(outFile,'    allocate(filterCX(%d,%d))\n',I-1,matrices.x.nTypes);
fprintf(outFile,'    allocate(filterRX(%d,%d,%d))\n\n',I,2*matrices.x.nRHSf-1,matrices.x.nTypes);

writeMatrix(outFile,'derivsAX',matrices.x.A);
writeMatrix(outFile,'derivsBX',matrices.x.B);
writeMatrix(outFile,'derivsCX',matrices.x.C);
writeMatrix(outFile,'derivsRX',matrices.x.R);

writeMatrix(outFile,'filterAX',matrices.x.Af);
writeMatrix(outFile,'filterBX',matrices.x.Bf);
writeMatrix(outFile,'filterCX',matrices.x.Cf);
writeMatrix(outFile,'filterRX',matrices.x.Rf);

fprintf(outFile,'    periodicX = %d\n\n',matrices.x.periodic);

if matrices.x.periodic
    fprintf(outFile,'    allocate(derivsDX(%d,%d))\n',I,matrices.x.nTypes);
    fprintf(outFile,'    allocate(filterDX(%d,%d))\n',I,matrices.x.nTypes);
    
    writeMatrix(outFile,'derivsDX',matrices.x.D);
    writeMatrix(outFile,'filterDX',matrices.x.Df);
end

%% For Y derivatives
fprintf(outFile,'    derivnRHSy = %d\n',matrices.y.nRHS);
fprintf(outFile,'    filternRHSy = %d\n',matrices.y.nRHSf);

fprintf(outFile,'    allocate(derivsAY(%d,%d))\n',J-1,matrices.y.nTypes);
fprintf(outFile,'    allocate(derivsBY(%d,%d))\n',J-matrices.y.periodic,matrices.y.nTypes);
fprintf(outFile,'    allocate(derivsCY(%d,%d))\n',J-1,matrices.y.nTypes);
fprintf(outFile,'    allocate(derivsRY(%d,%d,%d))\n\n',J,2*matrices.y.nRHS-1,matrices.y.nTypes);

fprintf(outFile,'    allocate(filterAY(%d,%d))\n',J-1,matrices.y.nTypes);
fprintf(outFile,'    allocate(filterBY(%d,%d))\n',J-matrices.y.periodic,matrices.y.nTypes);
fprintf(outFile,'    allocate(filterCY(%d,%d))\n',J-1,matrices.y.nTypes);
fprintf(outFile,'    allocate(filterRY(%d,%d,%d))\n\n',J,2*matrices.y.nRHSf-1,matrices.y.nTypes);

writeMatrix(outFile,'derivsAY',matrices.y.A);
writeMatrix(outFile,'derivsBY',matrices.y.B);
writeMatrix(outFile,'derivsCY',matrices.y.C);
writeMatrix(outFile,'derivsRY',matrices.y.R);

writeMatrix(outFile,'filterAY',matrices.y.Af);
writeMatrix(outFile,'filterBY',matrices.y.Bf);
writeMatrix(outFile,'filterCY',matrices.y.Cf);
writeMatrix(outFile,'filterRY',matrices.y.Rf);

fprintf(outFile,'    periodicY = %d\n\n',matrices.y.periodic);

if matrices.y.periodic
    fprintf(outFile,'    allocate(derivsDY(%d,%d))\n',J,matrices.y.nTypes);
    fprintf(outFile,'    allocate(filterDY(%d,%d))\n',J,matrices.y.nTypes);
    
    writeMatrix(outFile,'derivsDY',matrices.y.D);
    writeMatrix(outFile,'filterDY',matrices.y.Df);
end

%% For Z derivatives
if K > 1
    fprintf(outFile,'    derivnRHSz = %d\n',matrices.z.nRHS);
    fprintf(outFile,'    filternRHSz = %d\n',matrices.z.nRHSf);

    fprintf(outFile,'    allocate(derivsAZ(%d,%d))\n',K-1,matrices.z.nTypes);
    fprintf(outFile,'    allocate(derivsBZ(%d,%d))\n',K-matrices.z.periodic,matrices.z.nTypes);
    fprintf(outFile,'    allocate(derivsCZ(%d,%d))\n',K-1,matrices.z.nTypes);
    fprintf(outFile,'    allocate(derivsRZ(%d,%d,%d))\n\n',K,2*matrices.z.nRHS-1,matrices.z.nTypes);

    fprintf(outFile,'    allocate(filterAZ(%d,%d))\n',K-1,matrices.z.nTypes);
    fprintf(outFile,'    allocate(filterBZ(%d,%d))\n',K-matrices.z.periodic,matrices.z.nTypes);
    fprintf(outFile,'    allocate(filterCZ(%d,%d))\n',K-1,matrices.z.nTypes);
    fprintf(outFile,'    allocate(filterRZ(%d,%d,%d))\n\n',K,2*matrices.z.nRHSf-1,matrices.z.nTypes);

    writeMatrix(outFile,'derivsAZ',matrices.z.A);
    writeMatrix(outFile,'derivsBZ',matrices.z.B);
    writeMatrix(outFile,'derivsCZ',matrices.z.C);
    writeMatrix(outFile,'derivsRZ',matrices.z.R);

    writeMatrix(outFile,'filterAZ',matrices.z.Af);
    writeMatrix(outFile,'filterBZ',matrices.z.Bf);
    writeMatrix(outFile,'filterCZ',matrices.z.Cf);
    writeMatrix(outFile,'filterRZ',matrices.z.Rf);

    fprintf(outFile,'    periodicZ = %d\n\n',matrices.z.periodic);
    
    if matrices.z.periodic
        fprintf(outFile,'    allocate(derivsDZ(%d,%d))\n',K,matrices.z.nTypes);
        fprintf(outFile,'    allocate(filterDZ(%d,%d))\n',K,matrices.z.nTypes);

        writeMatrix(outFile,'derivsDZ',matrices.z.D);
        writeMatrix(outFile,'filterDZ',matrices.z.Df);
    end
    
end

%% Neumann Coefficients

fprintf(outFile,'\n    neumannLength = %d\n',length(matrices.neumannCoeffs));
fprintf(outFile,'    allocate(neumannCoeffs(%d))\n',length(matrices.neumannCoeffs));
fprintf(outFile,'    neumannCoeffs = (/');
for i = 1:length(matrices.neumannCoeffs)-1
    fprintf(outFile,'%.20fd0,',matrices.neumannCoeffs(i));
end
fprintf(outFile,'%.20fd0/)\n\n',matrices.neumannCoeffs(end));

fprintf(outFile,'\n    neumann2Length = %d\n',length(matrices.neumann2Coeffs));
fprintf(outFile,'    allocate(neumann2Coeffs(%d))\n',length(matrices.neumann2Coeffs));
fprintf(outFile,'    neumann2Coeffs = (/');
for i = 1:length(matrices.neumann2Coeffs)-1
    fprintf(outFile,'%.20fd0,',matrices.neumann2Coeffs(i));
end
fprintf(outFile,'%.20fd0/)\n',matrices.neumann2Coeffs(end));

% Tracked points
if ~isfield(mesh,'trackedPoints') || isempty(mesh.trackedPoints)
    fprintf(outFile,'    nTracked = 0\n');
else
    fprintf(outFile,'    nTracked = %d\n',size(mesh.trackedPoints,1));
    fprintf(outFile,'    allocate(indTracked(%d,3))\n',size(mesh.trackedPoints,1));
    
    indTracked = mesh.trackedPoints;
    
	% xTemp is a vector that only contains the nodes that were originally in the mesh in case extraRefinement is used
	% This is to make sure that the tracked point remain the same when using extraRefinement
	xTemp = nan(size(mesh.X));
	xTemp(1:(mesh.x.extraRefinement+1):end) = mesh.X(1:(mesh.x.extraRefinement+1):end);
	yTemp = nan(size(mesh.Y));
	yTemp(1:(mesh.y.extraRefinement+1):end) = mesh.Y(1:(mesh.y.extraRefinement+1):end);
	zTemp = nan(size(mesh.Z));
	zTemp(1:(mesh.z.extraRefinement+1):end) = mesh.Z(1:(mesh.z.extraRefinement+1):end);
	
    for i = 1:size(mesh.trackedPoints,1)
        [~,indTracked(i,1)] = min(abs(indTracked(i,1)-xTemp));
        [~,indTracked(i,2)] = min(abs(indTracked(i,2)-yTemp));
        [~,indTracked(i,3)] = min(abs(indTracked(i,3)-zTemp));
    end
	clear xTemp yTemp zTemp
    
    fprintf(outFile,'    indTracked = reshape((/');
    for n = 1:length(indTracked(:))-1
        fprintf(outFile, '%d,', indTracked(n));
    end
    fprintf(outFile, '%d', indTracked(end));
    fprintf(outFile,'/),shape(indTracked))\n');
end

fclose(outFile);
end

function writeMatrix(outFile,name,var)
    fprintf(outFile,'    %s = reshape((/',name);
    for n = 1:length(var(:))-1
        fprintf(outFile, '%.20fd0,', var(n));
    end
    fprintf(outFile, '%.20fd0', var(end));
    fprintf(outFile,'/),shape(%s))\n',name);
end
