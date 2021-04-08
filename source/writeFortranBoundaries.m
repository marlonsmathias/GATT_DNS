function writeFortranBoundaries(caseName,bi)

% This function writes the boundaryInfo.F90 file, it has the data for all boundary conditions for each domain slice

nProcs = length(bi);

outFile = fopen([caseName '/bin/boundaryInfo.F90'],'w');

fprintf(outFile, '    select case (nrank)\n');

vars = {'U','V','W','E','P'};

for i = 1:nProcs
    
    fprintf(outFile, '        case (%d)\n',i-1);
    
    % Dirichlet
    for j = 1:length(vars)
        switch j
            case 1
                n = bi{i}.nUd;
                ind = bi{i}.iUd;
                val = bi{i}.vUd;
            case 2
                n = bi{i}.nVd;
                ind = bi{i}.iVd;
                val = bi{i}.vVd;
            case 3
                n = bi{i}.nWd;
                ind = bi{i}.iWd;
                val = bi{i}.vWd;
            case 4
                n = bi{i}.nEd;
                ind = bi{i}.iEd;
                val = bi{i}.vEd;
            case 5
                n = bi{i}.nPd;
                ind = bi{i}.iPd;
                val = bi{i}.vPd;
        end
        fprintf(outFile, '            n%sd = %d\n',vars{j},n); % nUd = 1
        
        if n > 0
            fprintf(outFile, '            allocate(i%sd(%d,6))\n',vars{j},n); % allocate(iUd(1,6))
            fprintf(outFile, '            allocate(v%sd(%d))\n',vars{j},n); % allocate(vUd(1))

            fprintf(outFile, '            i%sd = reshape((/',vars{j}); %iUd = reshape((/
            for k = 1:6*n-1
                fprintf(outFile,'%d,',ind(k)); % 1,2,3,4,
            end
            fprintf(outFile,'%d/),shape(i%sd))\n',ind(end),vars{j}); % 5/),size(iUd))
            
            fprintf(outFile, '            v%sd = (/',vars{j}); %vUd = (/
            for k = 1:n-1
                fprintf(outFile,'%.20fd0,',val(k)); % 1.0,2.0,3.0,4.0,
            end
            fprintf(outFile,'%.20fd0/)\n\n',val(end)); % 5.0/)
            
        end
    end
    
    % Neumann
    for j = 1:length(vars)
        switch j
            case 1
                n = bi{i}.nUn;
                ind = bi{i}.iUn;
                dir = bi{i}.dUn;
            case 2
                n = bi{i}.nVn;
                ind = bi{i}.iVn;
                dir = bi{i}.dVn;
            case 3
                n = bi{i}.nWn;
                ind = bi{i}.iWn;
                dir = bi{i}.dWn;
            case 4
                n = bi{i}.nEn;
                ind = bi{i}.iEn;
                dir = bi{i}.dEn;
            case 5
                n = bi{i}.nPn;
                ind = bi{i}.iPn;
                dir = bi{i}.dPn;
        end
        fprintf(outFile, '            n%sn = %d\n',vars{j},n); % nUn = 1
        
        if n > 0
            fprintf(outFile, '            allocate(i%sn(%d,6))\n',vars{j},n); % allocate(iUn(1,6))
            fprintf(outFile, '            allocate(d%sn(%d))\n',vars{j},n); % allocate(dUn(1))

            fprintf(outFile, '            i%sn = reshape((/',vars{j}); %iUd = reshape((/
            for k = 1:6*n-1
                fprintf(outFile,'%d,',ind(k)); % 1,2,3,4,
            end
            fprintf(outFile,'%d/),shape(i%sn))\n',ind(end),vars{j}); % 5/),size(iUn))
            
            fprintf(outFile, '            d%sn = (/',vars{j}); %dUn = (/
            for k = 1:n-1
                fprintf(outFile,'%d,',dir(k)); % 1,2,3,4,
            end
            fprintf(outFile,'%d/)\n\n',dir(end)); % 5/)
            
        else
            fprintf(outFile, '            allocate(i%sn(1,6))\n',vars{j}); % allocate(iUn(1,6))
            fprintf(outFile, '            allocate(d%sn(1))\n\n',vars{j}); % allocate(dUn(1))
        end
    end
    
    % Second derivative
    for j = 1:length(vars)
        switch j
            case 1
                n = bi{i}.nUs;
                ind = bi{i}.iUs;
                dir = bi{i}.dUs;
            case 2
                n = bi{i}.nVs;
                ind = bi{i}.iVs;
                dir = bi{i}.dVs;
            case 3
                n = bi{i}.nWs;
                ind = bi{i}.iWs;
                dir = bi{i}.dWs;
            case 4
                n = bi{i}.nEs;
                ind = bi{i}.iEs;
                dir = bi{i}.dEs;
            case 5
                n = bi{i}.nPs;
                ind = bi{i}.iPs;
                dir = bi{i}.dPs;
        end
        fprintf(outFile, '            n%ss = %d\n',vars{j},n); % nUs = 1
        
        if n > 0
            fprintf(outFile, '            allocate(i%ss(%d,6))\n',vars{j},n); % allocate(iUs(1,6))
            fprintf(outFile, '            allocate(d%ss(%d))\n',vars{j},n); % allocate(dUs(1))

            fprintf(outFile, '            i%ss = reshape((/',vars{j}); %iUs = reshape((/
            for k = 1:6*n-1
                fprintf(outFile,'%d,',ind(k)); % 1,2,3,4,
            end
            fprintf(outFile,'%d/),shape(i%ss))\n',ind(end),vars{j}); % 5/),size(iUs))
            
            fprintf(outFile, '            d%ss = (/',vars{j}); %dUs = (/
            for k = 1:n-1
                fprintf(outFile,'%d,',dir(k)); % 1,2,3,4,
            end
            fprintf(outFile,'%d/)\n\n',dir(end)); % 5/)
            
        else
            fprintf(outFile, '            allocate(i%ss(1,6))\n',vars{j}); % allocate(iUs(1,6))
            fprintf(outFile, '            allocate(d%ss(1))\n\n',vars{j}); % allocate(dUs(1))
        end
    end
    
    % Corners
    fprintf(outFile, '            cN = %d\n',bi{i}.cN); % cN = 1
    fprintf(outFile, '            allocate(cL(%d,6))\n',max(bi{i}.cN,1)); % allocate(cL(1,6))
    fprintf(outFile, '            allocate(cD(%d,3))\n',max(bi{i}.cN,1)); % allocate(cD(1,3))
    fprintf(outFile, '            allocate(cAdiabatic(%d))\n',max(bi{i}.cN,1)); % allocate(cAdiabatic(1))
    
    if bi{i}.cN > 0
        fprintf(outFile, '            cL = reshape((/'); % cL = reshape((/
        for j = 1:6*bi{i}.cN-1
            fprintf(outFile,'%d,',bi{i}.cL(j)); % 1,2,3,4,
        end
        fprintf(outFile,'%d/),shape(cL))\n',bi{i}.cL(end)); % 5/),size(cL))
        
        fprintf(outFile, '            cD = reshape((/'); % cD = reshape((/
        for j = 1:3*bi{i}.cN-1
            fprintf(outFile,'%d,',bi{i}.cD(j)); % 1,2,3,4,
        end
        fprintf(outFile,'%d/),shape(cD))\n',bi{i}.cD(end)); % 5/),size(cD))
        
        fprintf(outFile, '            cAdiabatic = (/'); %cAdiabatic = (/
        for j = 1:bi{i}.cN-1
            fprintf(outFile,'%d,',bi{i}.adiabatic(j)); % 1.0,2.0,3.0,4.0,
        end
        fprintf(outFile,'%d/)\n\n',bi{i}.adiabatic(end)); % 5.0/)
    end
    
end

fprintf(outFile, '    end select\n');

fclose(outFile);

end
