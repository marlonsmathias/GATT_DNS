function disturbTypes = writeFortranDisturbances(caseName,bi,tridimensional)

% This function creates the runDisturbances.F90 and disturbances.F90 files
% runDisturbances.F90 is called when the boundary conditions are applied. It calls the appropriate routines for each domain slice
% runForcings.F90 is called at the end of the Navier Stokes equations routine. It calls the appropriate routines for each domain slice
% disturbances.F90 is a Fortran module file, which concatenates all the disturbances routines from the disturbances folder

nProcs = length(bi);

outFileDisturb = fopen([caseName '/bin/runDisturbances.F90'],'w');
if ~tridimensional
    outFileForcing = fopen([caseName '/bin/runForcings2D.F90'],'w');
    system(['touch ' caseName '/bin/runForcings3D.F90']);
else
    outFileForcing = fopen([caseName '/bin/runForcings3D.F90'],'w');
    system(['touch ' caseName '/bin/runForcings2D.F90']);
end

fprintf(outFileDisturb, '    select case (nrank)\n');
fprintf(outFileForcing, '    select case (nrank)\n');

disturbTypes = {};

for i = 1:nProcs
    fprintf(outFileDisturb, '        case (%d)\n',i-1);
    fprintf(outFileForcing, '        case (%d)\n',i-1);
    
    for j = 1:length(bi{i}.disturb)
        
        if bi{i}.disturb{j}.forcing
            outFile = outFileForcing;
        else
            outFile = outFileDisturb;
        end

        di = bi{i}.disturb{j};
        
        disturbTypes{end+1} = di.type; %#ok<AGROW>
        
        nx = di.ind(2)-di.ind(1)+1;
        ny = di.ind(4)-di.ind(3)+1;
        nz = di.ind(6)-di.ind(5)+1;
        
        fprintf(outFile,'            call %s(%d,%d,%d,(/',di.type,nx,ny,nz); % TS(nx,ny,nz,(/
        
        for k = 1:nx-1
            fprintf(outFile,'%.20fd0,',di.X(k)); %1.0,2.0,3.0,
        end
        fprintf(outFile,'%.20fd0/),(/',di.X(end)); %4.0/),(/
        for k = 1:ny-1
            fprintf(outFile,'%.20fd0,',di.Y(k)); %1.0,2.0,3.0,
        end
        fprintf(outFile,'%.20fd0/),(/',di.Y(end)); %4.0/),(/
        for k = 1:nz-1
            fprintf(outFile,'%.20fd0,',di.Z(k)); %1.0,2.0,3.0,
        end
        if bi{i}.disturb{j}.forcing
            fprintf(outFile,'%.20fd0/)',di.Z(end)); %4.0/)
        else
            fprintf(outFile,'%.20fd0/),t',di.Z(end)); %4.0/),t
        end
        
        if bi{i}.disturb{j}.forcing
            for k = 1:length(di.var)
                fprintf(outFile,',%s(%d:%d,%d:%d,%d:%d)',di.var(k),di.ind(1),di.ind(2),di.ind(3),di.ind(4),di.ind(5),di.ind(6)); %,U(1:2,3:4,5:6)
            end
            for k = 1:length(di.var)
                fprintf(outFile,',d%s(%d:%d,%d:%d,%d:%d)',di.var(k),di.ind(1),di.ind(2),di.ind(3),di.ind(4),di.ind(5),di.ind(6)); %,dU(1:2,3:4,5:6)
            end
        else
            for k = 1:length(di.var)
                fprintf(outFile,',%s(%d:%d,%d:%d,%d:%d)',di.var(k),di.ind(1),di.ind(2),di.ind(3),di.ind(4),di.ind(5),di.ind(6)); %,U(1:2,3:4,5:6)
            end
        end
        
        for k = 1:length(di.par) % Add extra parameters: loop each cell
			if isnumeric(di.par{k})
				if length(di.par{k}) == 1 % If cell is a single number
					fprintf(outFile,',%.20fd0',di.par{k}); % 1.0
				else % If cell is a vector
					fprintf(outFile,',(/'); %,(/
					for kk = 1:length(di.par{k})-1
						fprintf(outFile,'%.20fd0,',di.par{k}(kk)); %1.0,2.0,3.0,
					end
					fprintf(outFile,'%.20fd0/)',di.par{k}(end)); %4.0/)
				end
			else % If cell is a string
				fprintf(outFile,',''%s''',di.par{k}); % 1.0
			end
        end
        fprintf(outFile,')\n'); % )
        
    end
    
end

fprintf(outFileForcing, '    end select\n');
fprintf(outFileDisturb, '    end select\n');

fclose(outFileForcing);
fclose(outFileDisturb);

disturbTypes = unique(disturbTypes);

outFile = fopen([caseName '/bin/disturbances.F90'],'w');
fprintf(outFile,'    module disturbances\n\n    contains\n\n');

for i = 1:length(disturbTypes)
    if ~strcmp(disturbTypes{i},'holdInlet')
        sourceFile = fopen(['source/disturbances/' disturbTypes{i} '.F90'],'r');
    else
        sourceFile = fopen(['source/Fortran/' disturbTypes{i} '.F90'],'r');
    end
    line = fgetl(sourceFile);
    while ischar(line)
        fprintf(outFile,'%s\n',line);
        line = fgetl(sourceFile);
    end
    fclose(sourceFile);
    
    fprintf(outFile,'\n');
    
end

fprintf(outFile,'\n    end module\n');
fclose(outFile);

end
