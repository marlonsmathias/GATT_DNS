function runBatch(parName,stage)

if nargin == 1
	stage = [1 2 3];
end

if iscell(parName);
	for i = 1:length(parName)
		runBatch(parName{i},stage)
	end
	return
end

eval(parName)

maxNumCompThreads(p_row*p_col);

if any(stage == 1)
	[~,info] = runDNS(parName)
end

if any(stage == 2)
	if exist('info','var')
		mesh = info.mesh;
	else
		mesh = load([caseName '/mesh.mat'])
	end

	Ldomain = mesh.X(end) - mesh.X(1);
	%Twave = 2*4*Ldomain*Ma;
	c = 1/Ma;
	v1 = c+1;
	v2 = c-1;
	T1 = Ldomain/v1;
	T2 = Ldomain/v2;
	Twave = 4*(T1+T2);
	eta = 0.7;

	filesForMean = round(eta*Twave/time.qtimes); % Number of files for each mean
	spacingForMean = round(1.1*filesForMean); % Number of saved files between each mean
	nMeans = 10; % Number of loops
	tol = 1e-12;

	runDNSautomean_internal(parName,caseName,time.qtimes,time.control,filesForMean,spacingForMean,nMeans,tol)
end

if any(stage == 3)
	calcInstability_internal(caseName)
end


end

function calcInstability_internal(caseFolder)

%% Define parameters

% Case name
%caseFolder = 'Red600-Ma05-L10-D5-ppp13-refCav2h';
comment = 'SFD01';

% Method parameters
M = 2000; % Size of Krylov space
epsilon = 1e-6; % Disturbance magnitude
dT = 2;
singleBeta = true;
firstOrder = true;

removeBufferZone = false;

% SFD parameters (optional)
resumeSFD = true;

SFD.type = 2; % 0 = off, 1 = whole domain, 2 = buffer zone only;
SFD.Delta = inf;
SFD.X = 0.1;
SFD.applyY = true;

%SFD.applyY = true;

SFD.extraRegion = [];

% Domain parameters (optional)
%nz = 16;
%beta = pi;

% Domain decomposition (optional)
%global_p_row = 4;
%global_p_col = 1;

% Output parameters
nModes = 50; % Number of modes to be saved
removeCC = true;
nKrilov = M; % When the results files will be generated
saveEvery = 100; % When intermediary files will be saved for resuming

printEigsEvery = 10; % When intermediary eigenvalues will be printed
printEigsN = 6; % Number of eigenvalues printed

displayAllOutputs = false;

% Set folders for libraries
matlabDir = ''; % Leave empty for automatic directory
decompDir = '/usr/local/2decomp_fft';

% Parameters required by the DNS
logAll = false;
displayCompiling = false;
optimizeCode = true;
debugger = false;
profiler = false;
runningLST = true;

if displayAllOutputs
    suppressOutput = '';
else
    suppressOutput = ' > /dev/null 2>&1';
end

%% Read parameters and change some values
caseFile = 'parameters';
fprintf('Reading parameters file: ')
if exist([caseFolder '/' caseFile '.m'],'file') % First check case dir
    disp([caseFolder '/' caseFile '.m']);
    cd(caseFolder)
    eval(caseFile)
    cd ..
elseif exist([caseFolder '/bin/' caseFile '.m'],'file') % Then check bin dir
    disp([caseFolder '/bin/' caseFile '.m'])
    cd(caseFolder)
    cd bin
    eval(caseFile)
    cd ..
    cd ..
else % Now check for any .m in case dir
    allFiles = dir(caseFolder);
    parFiles = {};

    for i = 1:length(allFiles)
        name = allFiles(i).name;
        if endsWith(name,'.m')
            parFiles{end+1} = name; %#ok<AGROW>
        end
    end
    
    if length(parFiles) == 1
        disp([caseFolder '/' parFiles{1}]);
        cd(caseFolder)
        eval(parFiles{1}(1:end-2))
        cd ..
        
        
    else % Finally, check for any .m in bin dir
    
        allFiles = dir([caseFolder '/bin']);
        parFiles = {};

        for i = 1:length(allFiles)
            name = allFiles(i).name;
            if endsWith(name,'.m')
                parFiles{end+1} = name; %#ok<AGROW>
            end
        end
        
        if length(parFiles) == 1
            disp([caseFolder '/bin/' parFiles{1}]);
            cd(caseFolder)
            cd bin
            eval(parFiles{1}(1:end-2))
            cd ..
            cd ..
        else
        
            error('Parameters file not found')
        end
        
    end
    
end

caseName = caseFolder;
logAll = false;

if exist('beta','var') && beta > 0
    domain.zi = 0;
    domain.zf = 2*pi/beta;
    mesh.z.fixPeriodicDomainSize = true;
	if singleBeta
		mesh.z.n = 4;
	end
end

if exist('nz','var')
    mesh.z.n = nz;
end

if exist('SFD','var')
	extraVars = unpackStruct(SFD);
    for i = 1:length(extraVars)
        eval(['numMethods.SFD.' extraVars{i} ' = SFD.' extraVars{i} ';']);
    end
end


if exist('global_p_row','var')
    p_row = global_p_row;
end
if exist('global_p_col','var')
    p_col = global_p_col;
end

time.control = 'cfl';
time.qtimes = dT;
time.tmax = dT;

%% Create case name
if mesh.z.n == 1
    singleBeta = false;
end

caseNameInstability = [caseFolder '-dT' num2str(dT) '-eps' num2str(epsilon)];
if ~firstOrder
    caseNameInstability = [caseNameInstability '-SO'];
end
if mesh.z.n > 1 && ~singleBeta
    caseNameInstability = [caseNameInstability '-MB'];
end
if exist('nz','var')
    caseNameInstability = [caseNameInstability '-nz' num2str(nz)];
end
if exist('beta','var')
    caseNameInstability = [caseNameInstability '-beta' num2str(beta/pi) 'pi'];
end
if exist('comment','var') && ~isempty(comment)
    caseNameInstability = [caseNameInstability '-' comment];
end
fprintf('Case name: %s\n',caseNameInstability)

%% Run preprocessing
fprintf('Running preprocessor\n')
addpath source
addpath source/boundaries

% Clear previous system calls cache
[~,~] = system('');

% Create a copy of the original bin folder and mesh file to be restored after compling
cleanupObj = onCleanup(@()(cleanUpFunction(caseFolder)));
system(['cp -r ' caseFolder '/bin ' caseFolder '/bin_backup >/dev/null 2>&1']);
system(['cp ' caseFolder '/mesh.mat ' caseFolder '/mesh.mat.backup >/dev/null 2>&1']);
    
preprocessing

% Set maximum number of threads for Matlab processes
maxNumCompThreads(p_row*p_col);

% Fix nSave in parameters.F90 and set to 0
system(['sed -i ''s/    integer :: nSave = .*/    integer :: nSave = 0/'' ' caseFolder '/bin/parameters.F90']);

% Fix resumeSFD in parameters.F90
if resumeSFD
	system(['sed -i ''s/    integer :: resumeMeanFlow = .*/    integer :: resumeMeanFlow = 1/'' ' caseFolder '/bin/parameters.F90']);
end

%% Compile code
fprintf('Compiling code\n')
compileFortran

%% Create new folders and copy files
fprintf('Creating new folders and files\n')

if ~exist([caseFolder '/Instability'],'dir')
	instFolderName = ['/dev/shm/Instability' num2str(randi(1e8))];
    mkdir(instFolderName);
	system(['ln -s ' instFolderName ' ' caseFolder '/Instability']);
end

if ~exist([caseFolder '/Instability/' caseNameInstability],'dir')
    mkdir([caseFolder '/Instability/' caseNameInstability]);
else
    system(['rm ' caseFolder '/Instability/' caseNameInstability '/* -r']);
end

system(['cp ' caseFolder '/bin ' caseFolder '/Instability/' caseNameInstability '/bin -r']);

logFile = fopen([caseFolder '/Instability/' caseNameInstability '/log.txt'],'w');
fclose(logFile);

% Restore the original files
clear cleanupObj runningSFD

% Create new cleanUpObj to remove the Instability folder
cleanupObj = onCleanup(@()(cleanUpFunction(caseFolder,caseNameInstability)));

%% Read the base flow
fprintf('Reading base flow: ')

if exist([caseFolder '/baseflow.mat'],'file')
    baseflow = load([caseFolder '/baseflow.mat']);
	fprintf('%s/baseflow.mat\n',caseFolder)
else
    nStep = checkPreviousRun(caseFolder);
    baseflow = load(sprintf('%s/flow_%.10d.mat',caseFolder,nStep));
	fprintf('%s/flow_%.10d.mat\n',caseFolder,nStep);
end

% If needed, replicate flow to make it 3D
if mesh.z.n > 1 && size(baseflow.U,3) == 1
    baseflow.U = repmat(baseflow.U,[1 1 mesh.z.n]);
    baseflow.V = repmat(baseflow.V,[1 1 mesh.z.n]);
    baseflow.W = repmat(baseflow.W,[1 1 mesh.z.n]);
    baseflow.R = repmat(baseflow.R,[1 1 mesh.z.n]);
    baseflow.E = repmat(baseflow.E,[1 1 mesh.z.n]);
end

% Get the physical flow region
flowRegion = ~isnan(baseflow.U); % Get region outside of walls

if removeBufferZone
	flowRegion([1:mesh.x.buffer.i.n end-mesh.x.buffer.f.n+1:end],:,:) = false; % Remove buffer zones
	flowRegion(:,[1:mesh.y.buffer.i.n end-mesh.y.buffer.f.n+1:end],:) = false;
	flowRegion(:,:,[1:mesh.z.buffer.i.n end-mesh.z.buffer.f.n+1:end]) = false;
end

if mesh.z.n == 1
    N = 4*sum(flowRegion(:));
elseif ~singleBeta
    N = 5*sum(flowRegion(:));
else
    N = 5*sum(sum(flowRegion(:,:,1)));
end

%% Allocate variables
fprintf('Allocating variables\n');
zeta = zeros(N,M);
H = zeros(M,M);

%% Check for previous results and reload data
if exist([caseFolder '/results-' caseNameInstability '.mat'],'file')
    fprintf('Previous results file found, resuming\n')
    previousResults = matfile([caseFolder '/results-' caseNameInstability '.mat']);
    
    M0 = min(previousResults.M,M);
    zeta(:,1:M0) = previousResults.zeta(:,1:M0);
    H(1:M0,1:M0-1) = previousResults.H(1:M0,1:M0-1);
else
    M0 = 1;
    
    zeta(:,1) = calcInitialDisturb(mesh.X,mesh.Y,mesh.Z,singleBeta,flowRegion);
    
end

%% Set the initial SFD state
if resumeSFD
    t = 0;
    U = baseflow.U;
    V = baseflow.V;
    W = baseflow.W;
    R = baseflow.R;
    E = baseflow.E;
    save([caseFolder '/Instability/' caseNameInstability '/meanflowSFD.mat'],'U','V','W','R','E','t','-v7.3')
end

%% Iterate

if firstOrder % If using a first order approximation, compute the baseflow drift
	disp('Iteration 0')
    t = 0;
	
    U = baseflow.U;
    V = baseflow.V;
    W = baseflow.W;
    R = baseflow.R;
    E = baseflow.E;
    save([caseFolder '/Instability/' caseNameInstability '/flow_0000000000.mat'],'U','V','W','R','E','t','-v7.3')
    
    system(['(cd ' caseFolder '/Instability/' caseNameInstability '/bin && mpirun -np ' num2str(p_row*p_col) ' main ' caseNameInstability ') ' suppressOutput]);
    flowMinus = load([caseFolder '/Instability/' caseNameInstability '/flow_0000000001.mat']);
    delete([caseFolder '/Instability/' caseNameInstability '/flow_0000000001.mat'])
	
	UM = flow2vec(flowMinus,flowRegion,singleBeta);
end



tic
elapsedTime = 0;
for k = M0:M
    fprintf('\n');
    disp(['Iteration ' num2str(k) ' of ' num2str(M)]);
    
    % Get current disturbance
    disturbance = vec2flow(zeta(:,k)*epsilon*sqrt(N),flowRegion,singleBeta);
    
    % Save flows, run DNS and read results
    t = 0;
	
    U = baseflow.U + disturbance.U;
    V = baseflow.V + disturbance.V;
    W = baseflow.W + disturbance.W;
    R = baseflow.R + disturbance.R;
    E = baseflow.E + disturbance.E;
    save([caseFolder '/Instability/' caseNameInstability '/flow_0000000000.mat'],'U','V','W','R','E','t','-v7.3')
    
    system(['(cd ' caseFolder '/Instability/' caseNameInstability '/bin && mpirun -np ' num2str(p_row*p_col) ' main ' caseNameInstability ') ' suppressOutput]);
    flowPlus = load([caseFolder '/Instability/' caseNameInstability '/flow_0000000001.mat']);
    delete([caseFolder '/Instability/' caseNameInstability '/flow_0000000001.mat'])
    
	UP = flow2vec(flowPlus,flowRegion,singleBeta);
	
	if ~firstOrder
		U = baseflow.U - disturbance.U;
		V = baseflow.V - disturbance.V;
		W = baseflow.W - disturbance.W;
		R = baseflow.R - disturbance.R;
		E = baseflow.E - disturbance.E;
		save([caseFolder '/Instability/' caseNameInstability '/flow_0000000000.mat'],'U','V','W','R','E','t','-v7.3')
		
		system(['(cd ' caseFolder '/Instability/' caseNameInstability '/bin && mpirun -np ' num2str(p_row*p_col) ' main ' caseNameInstability ') ' suppressOutput]);
		flowMinus = load([caseFolder '/Instability/' caseNameInstability '/flow_0000000001.mat']);
		delete([caseFolder '/Instability/' caseNameInstability '/flow_0000000001.mat'])

		UM = flow2vec(flowMinus,flowRegion,singleBeta);

		% Compute matrices
		B = (UP - UM)/(2*epsilon*sqrt(N));
	else
		B = (UP - UM)/(epsilon*sqrt(N));
	end
    
    uPrime = B;

    for j = 1:k
        H(j,k) = zeta(:,j)'*B;
        uPrime = uPrime - H(j,k)*zeta(:,j);
    end
    
    if k<M
        H(k+1,k) = norm(uPrime);
        zeta(:,k+1) = uPrime/H(k+1,k);
    end
    
    if any(isnan(H(:))) || any(isinf(H(:)))
        error('DNS solution has failed')
    end
    
    %% Save results to file
    if mod(k,saveEvery) == 0 && k ~= M0
        actualM = M;
		actualH = H;
		actualZeta = zeta;
		
        M = k;
		
		H(M+1:end,:) = [];
		H(:,M+1:end) = [];
		zeta(:,M+1:end) = [];
		
        disp(['Saving results to results-' caseNameInstability '.mat'])
        save([caseFolder '/results-' caseNameInstability '.mat'],'zeta','H','M','-v7.3')
		
		zeta = actualZeta;
		H = actualH;
        M = actualM;
		
		clear actualZeta actualM actualM
    end

    %% Compute modes
    if any(k == nKrilov)
        disp('Computing modes')

        [psiH, lambdaH] = eig(H(1:k,1:k));
        lambdaB = diag(lambdaH);

        lambdaA = log(lambdaB)/dT;
        [~,index] = sort(real(lambdaA),'descend');
        lambdaA = lambdaA(index);

        psiH = psiH(:,index);

        psiA = zeta(:,1:k)*psiH(:,1:min(nModes,k));
        
        modes = vec2modes(psiA,flowRegion,singleBeta,min(nModes,k),lambdaA,removeCC);
        modes.lambda = lambdaA;
        modes.X = mesh.X;
        modes.Y = mesh.Y;
        modes.Z = mesh.Z;
        
        %% Save modes
        disp(['Saving modes to modes-' caseNameInstability '-M' num2str(k) '.mat'])
        save([caseFolder '/modes-' caseNameInstability '-M' num2str(k) '.mat'],'-struct','modes','-v7.3')
        
    end
    
    %% Print remaining time
    elapsedTime(end+1) = toc;
    remainingTime = (elapsedTime(end)-elapsedTime(1))/(length(elapsedTime)-1)*(M-k);
	if length(elapsedTime) > 100
		elapsedTime(1) = [];
	end
    disp(['Elapsed time: ' sec2str(elapsedTime(end))])
    disp(['Estimated time remaining: ' sec2str(remainingTime)])
	
	%% Print leading eigenvalues
	if mod(k,printEigsEvery) == 0
		leadingEigs = eigs(H,min(k,printEigsN),'largestabs','tol',1e-4,'MaxIterations',max(300,2*size(H,1)),'FailureTreatment','keep');
		leadingEigs = log(leadingEigs')/dT;
		disp(['Leading eigenvalues : ' num2str(leadingEigs)]);
		if exist('leadingEigsPrev','var')
			disp(['Change : ' num2str(abs((leadingEigs(1:length(leadingEigsPrev))-leadingEigsPrev)./leadingEigs(1:length(leadingEigsPrev))))]);
		end
		leadingEigsPrev = leadingEigs;
	end
    
end
end

%% Extra functions
function vec = flow2vec(flow,flowRegion,singleBeta)
    threeD = size(flowRegion,3) > 1;
    
    if ~threeD % 2D flow
    
        vec = [flow.U(flowRegion); flow.V(flowRegion); flow.R(flowRegion); flow.E(flowRegion)];
        
    elseif ~singleBeta % Fully 3D flow
    
        vec = [flow.U(flowRegion); flow.V(flowRegion); flow.W(flowRegion); flow.R(flowRegion); flow.E(flowRegion)];
    
    else % Single beta, 3D flow
        
        if size(flow.U,3) > 1
            nz = size(flowRegion,3);
            Ztemp = linspace(0,2*pi,nz+1);
            Ztemp(end) = [];
            Ztemp = permute(Ztemp,[1 3 2]);

            cosZ = cos(Ztemp)*2/nz;
            sinZ = sin(Ztemp)*2/nz;
            
            flow.U = sum(cosZ.*flow.U,3);
            flow.V = sum(cosZ.*flow.V,3);
            flow.W = sum(sinZ.*flow.W,3);
            flow.R = sum(cosZ.*flow.R,3);
            flow.E = sum(cosZ.*flow.E,3);
        end
        
        flowRegion = flowRegion(:,:,1);
        
        vec = [flow.U(flowRegion); flow.V(flowRegion); flow.W(flowRegion); flow.R(flowRegion); flow.E(flowRegion)];
    
    end
    
end

function flow = vec2flow(vec,flowRegion,singleBeta)
    
    [nx, ny, nz] = size(flowRegion);
    N = sum(flowRegion(:));
    
    threeD = nz > 1;
    
    if ~threeD % 2D flow
    
        flow.U = zeros(nx,ny);
        flow.V = zeros(nx,ny);
        flow.W = zeros(nx,ny);
        flow.R = zeros(nx,ny);
        flow.E = zeros(nx,ny);
        
        flow.U(flowRegion) = vec(1:N);
        flow.V(flowRegion) = vec(N+1:2*N);
        flow.R(flowRegion) = vec(2*N+1:3*N);
        flow.E(flowRegion) = vec(3*N+1:4*N);
        
    elseif ~singleBeta % Fully 3D flow
    
        flow.U = zeros(nx,ny,nz);
        flow.V = zeros(nx,ny,nz);
        flow.W = zeros(nx,ny,nz);
        flow.R = zeros(nx,ny,nz);
        flow.E = zeros(nx,ny,nz);
        
        flow.U(flowRegion) = vec(1:N);
        flow.V(flowRegion) = vec(N+1:2*N);
        flow.W(flowRegion) = vec(2*N+1:3*N);
        flow.R(flowRegion) = vec(3*N+1:4*N);
        flow.E(flowRegion) = vec(4*N+1:5*N);
    
    else % Single beta, 3D flow
        
        flow.U = zeros(nx,ny);
        flow.V = zeros(nx,ny);
        flow.W = zeros(nx,ny);
        flow.R = zeros(nx,ny);
        flow.E = zeros(nx,ny);
        
        N = round(N/nz);
        
        flow.U(flowRegion(:,:,1)) = vec(1:N);
        flow.V(flowRegion(:,:,1)) = vec(N+1:2*N);
        flow.W(flowRegion(:,:,1)) = vec(2*N+1:3*N);
        flow.R(flowRegion(:,:,1)) = vec(3*N+1:4*N);
        flow.E(flowRegion(:,:,1)) = vec(4*N+1:5*N);
    
        Ztemp = linspace(0,2*pi,nz+1);
        Ztemp(end) = [];
        Ztemp = permute(Ztemp,[1 3 2]);
        
        cosZ = cos(Ztemp);
        sinZ = sin(Ztemp);
        
        flow.U = flow.U.*cosZ;
        flow.V = flow.V.*cosZ;
        flow.W = flow.W.*sinZ;
        flow.R = flow.R.*cosZ;
        flow.E = flow.E.*cosZ;
    
    end
    
end

function modes = vec2modes(vec,flowRegion,singleBeta,nModes,lambda,removeCC)
    
    [nx, ny, nz] = size(flowRegion);
    N = sum(flowRegion(:));
    
    threeD = nz > 1;
    
    if ~threeD % 2D flow
    
        modes.U = nan(nx,ny,1,nModes);
        modes.V = nan(nx,ny,1,nModes);
        modes.W = nan(nx,ny,1,nModes);
        modes.R = nan(nx,ny,1,nModes);
        modes.E = nan(nx,ny,1,nModes);
        
        U = nan(nx,ny);
        V = nan(nx,ny);
        R = nan(nx,ny);
        E = nan(nx,ny);
        
        for i = 1:nModes
			if imag(lambda(i))>=0
				vec(:,i) = vec(:,i)/max(abs(vec(:,i)));
				
				U(flowRegion) = vec(1:N,i);
				V(flowRegion) = vec(N+1:2*N,i);
				R(flowRegion) = vec(2*N+1:3*N,i);
				E(flowRegion) = vec(3*N+1:4*N,i);
				
				modes.U(:,:,1,i) = U;
				modes.V(:,:,1,i) = V;
				modes.R(:,:,1,i) = R;
				modes.E(:,:,1,i) = E;
			end
        end
        
    elseif ~singleBeta % Fully 3D flow
    
        modes.U = nan(nx,ny,nz,nModes);
        modes.V = nan(nx,ny,nz,nModes);
        modes.W = nan(nx,ny,nz,nModes);
        modes.R = nan(nx,ny,nz,nModes);
        modes.E = nan(nx,ny,nz,nModes);
        
        U = nan(nx,ny,nz);
        V = nan(nx,ny,nz);
        W = nan(nx,ny,nz);
        R = nan(nx,ny,nz);
        E = nan(nx,ny,nz);
        
        for i = 1:nModes
            if imag(lambda(i))>=0
				vec(:,i) = vec(:,i)/max(abs(vec(:,i)));
				
				U(flowRegion) = vec(1:N,i);
				V(flowRegion) = vec(N+1:2*N,i);
				W(flowRegion) = vec(2*N+1:3*N,i);
				R(flowRegion) = vec(3*N+1:4*N,i);
				E(flowRegion) = vec(4*N+1:5*N,i);
				
				modes.U(:,:,:,i) = U;
				modes.V(:,:,:,i) = V;
				modes.W(:,:,:,i) = W;
				modes.R(:,:,:,i) = R;
				modes.E(:,:,:,i) = E;
			end
        end
    
    else % Single beta, 3D flow
        
        N = round(N/nz);
        flowRegion = flowRegion(:,:,1);
        
        modes.U = nan(nx,ny,1,nModes);
        modes.V = nan(nx,ny,1,nModes);
        modes.W = nan(nx,ny,1,nModes);
        modes.R = nan(nx,ny,1,nModes);
        modes.E = nan(nx,ny,1,nModes);
        
        U = nan(nx,ny);
        V = nan(nx,ny);
        W = nan(nx,ny);
        R = nan(nx,ny);
        E = nan(nx,ny);
        
        for i = 1:nModes
            if imag(lambda(i))>=0
				vec(:,i) = vec(:,i)/max(abs(vec(:,i)));
				
				U(flowRegion) = vec(1:N,i);
				V(flowRegion) = vec(N+1:2*N,i);
				W(flowRegion) = vec(2*N+1:3*N,i);
				R(flowRegion) = vec(3*N+1:4*N,i);
				E(flowRegion) = vec(4*N+1:5*N,i);
				
				modes.U(:,:,1,i) = U;
				modes.V(:,:,1,i) = V;
				modes.W(:,:,1,i) = W;
				modes.R(:,:,1,i) = R;
				modes.E(:,:,1,i) = E;
			end
        end
    
    end
    
end

function zeta = calcInitialDisturb(X,Y,Z,singleBeta,flowRegion)

    [~,x0i] = min(gradient(X));
    [~,y0i] = min(gradient(Y));
    [~,z0i] = min(gradient(Z));

    x0 = X(x0i);
    y0 = Y(y0i);
    z0 = Z(z0i);
    
    alphaX = 10/(X(end)-X(1)).^2;
    alphaY = 10/(Y(end)-Y(1)).^2;
    alphaZ = 10/(Z(end)-Z(1)).^2;
    
    if length(Z) == 1
        eta = exp(-(alphaX*(X'-x0).^2+alphaY*(Y-y0).^2));
        
        flow.U = eta;
        flow.V = eta;
        flow.W = eta;
        flow.R = eta;
        flow.E = eta;
        
        zeta = flow2vec(flow,flowRegion,singleBeta);
        
    elseif ~singleBeta
        eta = exp(-(alphaX*(X'-x0).^2+alphaY*(Y-y0).^2+alphaZ*(permute(Z,[1 3 2])-z0).^2));
        
        flow.U = eta;
        flow.V = eta;
        flow.W = eta;
        flow.R = eta;
        flow.E = eta;
        
        zeta = flow2vec(flow,flowRegion,singleBeta);
        
    else
        eta = exp(-(alphaX*(X'-x0).^2+alphaY*(Y-y0).^2));
        
        flow.U = eta;
        flow.V = eta;
        flow.W = eta;
        flow.R = eta;
        flow.E = eta;
        
        zeta = flow2vec(flow,flowRegion,singleBeta);
    end
    
    zeta = zeta/norm(zeta);
    
end

function timeString = sec2str(sec)
    sec = round(sec);
    if sec < 60
        timeString = sprintf('%ds',sec);
        return
    else
        min = floor(sec/60);
        sec = mod(sec,60);
    end
    if min < 60
        timeString = sprintf('%dmin %02ds',min,sec);
        return
    else
        hour = floor(min/60);
        min = mod(min,60);
    end
    if hour < 24
        timeString = sprintf('%dh %02dmin %02ds',hour,min,sec);
        return
    else
        day = floor(hour/24);
        hour = mod(hour,24);
    end
    timeString = sprintf('%dd %02dh %02dmin %02ds',day,hour,min,sec);
end

function varList = unpackStruct(structure)
    varList = {};
    varNames = fieldnames(structure);
    for i = 1:length(varNames)
        varName = varNames{i};
        if isstruct(structure.(varName))
            varList2 = unpackStruct(structure.(varName));
            for j = 1:length(varList2)
                varList{end+1} = [varName '.' varList2{j}];
            end
        else
            varList{end+1} = varName;
        end
    end
end

function cleanUpFunction(caseFolder,caseNameInstability)
	if nargin == 1 % Restore the original files
		system(['mv ' caseFolder '/mesh.mat.backup ' caseFolder '/mesh.mat >/dev/null 2>&1']);
		system(['rm -r ' caseFolder '/bin >/dev/null 2>&1']);
		system(['mv ' caseFolder '/bin_backup ' caseFolder '/bin >/dev/null 2>&1']);
	else % Clean up instability folder
		system(['rm -r ' caseFolder '/Instability/' caseNameInstability]);
		temp = dir([caseFolder '/Instability']);
		if length(temp) == 2 % Check if Instability dir is empty and erase it
			system(['rm -r "$(readlink -f ' caseFolder '/Instability )"']);
			system(['rm -r ' caseFolder '/Instability']);
		end
	end
end

function runDNSautomean_internal(parName,caseName,qtimes,timecontrol,filesForMean,spacingForMean,nMeans,tol)
%caseName = 'Red600-Ma05-L10-D5-ppp13-refCav2d';

%filesForMean = 139; % Number of files for each mean
%spacingForMean = 150; % Number of saved files between each mean
%nMeans = 3; % Number of loops

nMean = 1;

while true

	[~,changesStr] = system(['tail -1 ' caseName '/log.txt']);
	changes = str2num(changesStr);
	change = max(changes(6:10));
	
	if nMean > nMeans || change < tol
		break
	end
	
	fprintf('Computing means loop: %d\n Current change: %g\n', nMean, change)
	nMean = nMean + 1;

    caseFiles = {};
    allFiles = dir(caseName);
    for i = 1:length(allFiles)
        name = allFiles(i).name;
        if length(name) == 19 && ~isempty(regexp(name,'flow_\d*.mat','once')) % Check the file name
            caseFiles{end+1} = name;
        end
    end
    
    caseFiles = caseFiles(end-filesForMean+1:end);
    
    lastStep = sscanf(caseFiles{end},'flow_%d.mat');
    
    U = 0; V = 0; W = 0; R = 0; E = 0;
    for i = 1:filesForMean
        
        disp(caseFiles{i})
        current = load([caseName '/' caseFiles{i}]);

        U = U + current.U;
        V = V + current.V;
        W = W + current.W;
        R = R + current.R;
        E = E + current.E;
    end

    U = U/filesForMean;
    V = V/filesForMean;
    W = W/filesForMean;
    R = R/filesForMean;
    E = E/filesForMean;
    t = current.t;
    
    tstr = sprintf('%s/flow_%010d.mat',caseName,lastStep+1);
    save(tstr,'t','U','V','W','R','E')
    
    % Rename meanflowSFD
    system(['mv ' caseName '/meanflowSFD.mat ' caseName '/meanflowSFD_' num2str(lastStep) '.mat 2>/dev/null']);
    
    % Set new parameters
    switch timecontrol
        case 'cfl'
			extraPars.time.tmax = t + spacingForMean*qtimes;
        case 'dt'
			extraPars.time.tmax = lastStep + spacingForMean*qtimes;
    end
	
    % Write to log2 file
    logFile2 = fopen([caseName '/bin/log2.txt'],'a');
    fprintf(logFile2,'Mean flow calculated from %s to %s\n', caseFiles{1}, caseFiles{end});
    fclose(logFile2);
	
    % Call DNS
	runDNS(parName,extraPars)
		
end

end