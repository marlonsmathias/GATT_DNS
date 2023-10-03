% This script compiles the Fortran files

%% Create extra makefile with the directories that are specific to this run
outFile = fopen([caseName '/bin/makefile_extra'],'w');

if isempty(matlabDir)
    matlabDir = matlabroot;
end

fprintf(outFile,'MATROOT = %s\n',matlabDir);

fprintf(outFile,'DECOMPDIR = %s\n',decompDir);

if optimizeCode && ~debugger % Optimization options
    fprintf(outFile,'ARGS += -O5 -fcheck=all -fno-finite-math-only -march=native\n');
end

if debugger % Debugging options
    fprintf(outFile,'ARGS += -O0 -g -fbounds-check\n');
elseif profiler % Profiling options
    fprintf(outFile,'ARGS += -g -pg\n');
end

fclose(outFile);

%% Run make
if displayCompiling
	supressOutput = '';
else
	supressOutput = ' >/dev/null'; % This supresses the compiler output
end

if exist([caseName '/bin/main'],'file') % Remove the main binary to force recompiling
    delete([caseName '/bin/main'])
end

status = system(['cd ' caseName '/bin && make --makefile=../../source/Fortran/makefile' supressOutput]); % Call make

if status ~= 0
    error('Fortran compiling has failed')
end
