caseName = 'ReD3000-Ddelta5-Ma03-LDinf';

filesForMean = 90; % Number of files for each mean
spacingForMean = 100; % Number of saved files between each mean
nMeans = 20; % Number of loops
tol = 1e-12; % Tolerace for stopping

parFile = fopen([caseName '/bin/parameters.m']);

% Load base parameters file
parameters = {};
line = fgetl(parFile);
while ischar(line)

    if contains(line,'time.qtimes') || contains(line,'time.control') % get qtimes
        eval(line);
    end
        
	parameters{end+1} = line;
    line = fgetl(parFile);
end
fclose(parFile);

md = java.security.MessageDigest.getInstance('MD5');
hash = sprintf('%2.2x', typecast(md.digest(uint8(caseName)), 'uint8')');
parFileName = ['parameters_automean_' hash];

system(['touch ' parFileName '.m']);
cleanupObj = onCleanup(@()(cleanUpFunction(parFileName)));

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
    
    % Write new parameters file
    fPar = fopen([parFileName '.m'],'w');
    for i = 1:length(parameters)
        fprintf(fPar, '%s\n', parameters{i});
    end
    
    fprintf(fPar,'\ncaseName = ''%s'';\n', caseName);
    switch time.control
        case 'cfl'
            fprintf(fPar,'\ntime.tmax = %f;\n', t + spacingForMean*time.qtimes);
        case 'dt'
            fprintf(fPar,'\ntime.tmax = %s;\n', lastStep + spacingForMean*time.qtimes);
    end
    
    fclose(fPar);
	
    % Write to log2 file
    logFile2 = fopen([caseName '/bin/log2.txt'],'a');
    fprintf(logFile2,'Mean flow calculated from %s to %s\n', caseFiles{1}, caseFiles{end});
    fclose(logFile2);
	
    % Call DNS
	rehash PATH % This forces Matlab to reread the parameters file
	runDNS(parFileName)
		
end

clear cleanupObj

function cleanUpFunction(parFileName)
	delete([parFileName '.m'])
end
