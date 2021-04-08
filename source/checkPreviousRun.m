function [nStep, nx, ny, nz] = checkPreviousRun(caseName)

% This function checks for previous run files
% caseName is the folder to be checked
% nStep is the last time step found
% nx, ny, and nz are the size of the mesh in the saved file
% If no files are found, empty arrays are returned

allFiles = dir(caseName); % List all files

caseFiles = {}; % Flow file will be placed here

for i = 1:length(allFiles)
    name = allFiles(i).name;
    if length(name) == 19 && ~isempty(regexp(name,'flow_\d*.mat','once')) % Check the file name
        caseFiles{end+1} = name; %#ok<AGROW>
    end
end

if isempty(caseFiles)
    nStep = [];
    nx = [];
    ny = [];
    nz = [];
    return
end

nSteps = cellfun(@(name)(sscanf(name,'flow_%d.mat')),caseFiles);

nStep = max(nSteps);

if nargout > 1
	fileObject = matfile(sprintf('%s/flow_%.10d.mat',caseName,nStep));
	[nx, ny, nz] = size(fileObject.U);
end

end
