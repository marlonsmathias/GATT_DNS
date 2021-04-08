function bi = initBoundaries(boundary,mesh,domainSlicesY,domainSlicesZ,p_row,p_col)

% This fuction initializes the boundary conditions to the domain
% Before running the DNS, it should be called with the boundary structure
% as input

% Variable names are as follows:
% First letter - What it represents
%   n is the number of boundaries of that type
%   i is the index this boundary is located
%   v is the value a Dirichlet boundary will be set to
%   d is the direction a null Neumann boundary will be applied
% Second letter - The flow variable it applies to
% Third letter - The type of boundary: d for Dirichlet and n for Neumann

% One output is created for each domain slice
    
% Loop through all boundaries and organize them in the correct variables
biG.nUd=0; biG.nVd=0; biG.nWd=0; biG.nPd=0; biG.nEd=0;
biG.nUn=0; biG.nVn=0; biG.nWn=0; biG.nPn=0; biG.nEn=0;
biG.nUs=0; biG.nVs=0; biG.nWs=0; biG.nPs=0; biG.nEs=0;
biG.iUd = []; biG.iVd = []; biG.iWd = []; biG.iPd = []; biG.iEd = [];
biG.iUn = []; biG.iVn = []; biG.iWn = []; biG.iPn = []; biG.iEn = [];
biG.iUs = []; biG.iVs = []; biG.iWs = []; biG.iPs = []; biG.iEs = [];
biG.vUd = []; biG.vVd = []; biG.vWd = []; biG.vPd = []; biG.vEd = [];
biG.dUn = []; biG.dVn = []; biG.dWn = []; biG.dPn = []; biG.dEn = [];
biG.dUs = []; biG.dVs = []; biG.dWs = []; biG.dPs = []; biG.dEs = [];

directionOrder = {'xi','xf','yi','yf','zi','zf'};
for i = 1:length(boundary.val)
    switch boundary.type{i}
        case 'dir'
            switch boundary.var{i}
                case 'u'
                    biG.nUd = biG.nUd + 1;
                    biG.iUd(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)]; %#ok<*AGROW>
                    biG.vUd(end+1) = boundary.val(i);
                case 'v'
                    biG.nVd = biG.nVd + 1;
                    biG.iVd(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.vVd(end+1) = boundary.val(i);
                case 'w'
                    biG.nWd = biG.nWd + 1;
                    biG.iWd(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.vWd(end+1) = boundary.val(i);
                case 'p'
                    biG.nPd = biG.nPd + 1;
                    biG.iPd(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.vPd(end+1) = boundary.val(i);
                case 'e'
                    biG.nEd = biG.nEd + 1;
                    biG.iEd(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.vEd(end+1) = boundary.val(i);
            end
        case 'neu'
            switch boundary.var{i}
                case 'u'
                    biG.nUn = biG.nUn + 1;
                    biG.iUn(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.dUn(end+1) = find(strcmp(directionOrder,boundary.dir{i}));
                case 'v'
                    biG.nVn = biG.nVn + 1;
                    biG.iVn(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.dVn(end+1) = find(strcmp(directionOrder,boundary.dir{i}));
                case 'w'
                    biG.nWn = biG.nWn + 1;
                    biG.iWn(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.dWn(end+1) = find(strcmp(directionOrder,boundary.dir{i}));
                case 'p'
                    biG.nPn = biG.nPn + 1;
                    biG.iPn(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.dPn(end+1) = find(strcmp(directionOrder,boundary.dir{i}));
                case 'e'
                    biG.nEn = biG.nEn + 1;
                    biG.iEn(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.dEn(end+1) = find(strcmp(directionOrder,boundary.dir{i}));
            end
        case 'sec'
            switch boundary.var{i}
                case 'u'
                    biG.nUs = biG.nUs + 1;
                    biG.iUs(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.dUs(end+1) = find(strcmp(directionOrder,boundary.dir{i}));
                case 'v'
                    biG.nVs = biG.nVs + 1;
                    biG.iVs(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.dVs(end+1) = find(strcmp(directionOrder,boundary.dir{i}));
                case 'w'
                    biG.nWs = biG.nWs + 1;
                    biG.iWs(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.dWs(end+1) = find(strcmp(directionOrder,boundary.dir{i}));
                case 'p'
                    biG.nPs = biG.nPs + 1;
                    biG.iPs(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.dPs(end+1) = find(strcmp(directionOrder,boundary.dir{i}));
                case 'e'
                    biG.nEs = biG.nEs + 1;
                    biG.iEs(end+1,:) = [boundary.xi(i) boundary.xf(i) boundary.yi(i) boundary.yf(i) boundary.zi(i) boundary.zf(i)];
                    biG.dEs(end+1) = find(strcmp(directionOrder,boundary.dir{i}));
            end
    end
end

% Get information on corners
biG.cL = boundary.corners.limits;
biG.cD = boundary.corners.dir;
biG.adiabatic = boundary.corners.adiabatic;
biG.cN = size(biG.cL,1);

neumannLength = length(boundary.neumannCoeffs);
neumann2Length = length(boundary.neumann2Coeffs);

% Set gamma-1 value for converting pressure to density
biG.gamma1 = boundary.gamma - 1;

% Find regions inside walls
%biG.wR = ~boundary.flowRegion;
biG.E0 = boundary.E0;

%% Split boundaries for different processors

bi = cell(1,p_row*p_col);
for j = 1:p_row
    for k = 1:p_col
        nProc = k + (j-1)*p_col;
        
        biL = biG;
        
        Ji = domainSlicesY(1,j);
        Jf = domainSlicesY(2,j);
        Ki = domainSlicesZ(1,k);
        Kf = domainSlicesZ(2,k);
        
        [biL.iUd, biL.nUd, biL.vUd] = limitIndices(biL.iUd, biL.nUd, biL.vUd,'d',Ji,Jf,Ki,Kf,0);
        [biL.iVd, biL.nVd, biL.vVd] = limitIndices(biL.iVd, biL.nVd, biL.vVd,'d',Ji,Jf,Ki,Kf,0);
        [biL.iWd, biL.nWd, biL.vWd] = limitIndices(biL.iWd, biL.nWd, biL.vWd,'d',Ji,Jf,Ki,Kf,0);
        [biL.iPd, biL.nPd, biL.vPd] = limitIndices(biL.iPd, biL.nPd, biL.vPd,'d',Ji,Jf,Ki,Kf,0);
        [biL.iEd, biL.nEd, biL.vEd] = limitIndices(biL.iEd, biL.nEd, biL.vEd,'d',Ji,Jf,Ki,Kf,0);
        
        [biL.iUn, biL.nUn, biL.dUn] = limitIndices(biL.iUn, biL.nUn, biL.dUn,'n',Ji,Jf,Ki,Kf,neumannLength);
        [biL.iVn, biL.nVn, biL.dVn] = limitIndices(biL.iVn, biL.nVn, biL.dVn,'n',Ji,Jf,Ki,Kf,neumannLength);
        [biL.iWn, biL.nWn, biL.dWn] = limitIndices(biL.iWn, biL.nWn, biL.dWn,'n',Ji,Jf,Ki,Kf,neumannLength);
        [biL.iPn, biL.nPn, biL.dPn] = limitIndices(biL.iPn, biL.nPn, biL.dPn,'n',Ji,Jf,Ki,Kf,neumannLength);
        [biL.iEn, biL.nEn, biL.dEn] = limitIndices(biL.iEn, biL.nEn, biL.dEn,'n',Ji,Jf,Ki,Kf,neumannLength);
        
        [biL.iUs, biL.nUs, biL.dUs] = limitIndices(biL.iUs, biL.nUs, biL.dUs,'n',Ji,Jf,Ki,Kf,neumann2Length);
        [biL.iVs, biL.nVs, biL.dVs] = limitIndices(biL.iVs, biL.nVs, biL.dVs,'n',Ji,Jf,Ki,Kf,neumann2Length);
        [biL.iWs, biL.nWs, biL.dWs] = limitIndices(biL.iWs, biL.nWs, biL.dWs,'n',Ji,Jf,Ki,Kf,neumann2Length);
        [biL.iPs, biL.nPs, biL.dPs] = limitIndices(biL.iPs, biL.nPs, biL.dPs,'n',Ji,Jf,Ki,Kf,neumann2Length);
        [biL.iEs, biL.nEs, biL.dEs] = limitIndices(biL.iEs, biL.nEs, biL.dEs,'n',Ji,Jf,Ki,Kf,neumann2Length);

        values = [biL.cD biL.adiabatic];
        [biL.cL, biL.cN, values] = limitIndices(biL.cL, biL.cN, values,'c',Ji,Jf,Ki,Kf,neumannLength);
        biL.cD = values(:,1:3);
        biL.adiabatic = values(:,4);
        
        bi{nProc} = biL;
    end
end

%% Split disturbances across processors
for j = 1:p_row
    for k = 1:p_col
        nProc = k + (j-1)*p_col;
        disturb = {};
        for i = 1:length(boundary.disturb)
            if ~isempty(boundary.disturb{i})
                ind = boundary.disturb{i}.ind;
                ind(3:6) = [max(ind(3),domainSlicesY(1,j)) min(ind(4),domainSlicesY(2,j)) max(ind(5),domainSlicesZ(1,k)) min(ind(6),domainSlicesZ(2,k))];

                if ind(:,3) <= ind(:,4) && ind(:,5) <= ind(:,6)
                    disturb{end+1} = boundary.disturb{i};
                    disturb{end}.ind = ind;
                    disturb{end}.X = mesh.X(ind(1):ind(2));
                    disturb{end}.Y = mesh.Y(ind(3):ind(4));
                    disturb{end}.Z = mesh.Z(ind(5):ind(6));
                end
            end
        end
        
        bi{nProc}.disturb = disturb;
        
    end
end

end

function [ind, n, vd] = limitIndices(ind,n,vd,type,Ji,Jf,Ki,Kf,neumannLength)
    if n == 0
        return
    end
    for i = 1:n
        ind(i,[3 4 5 6]) = [max(ind(i,3),Ji) min(ind(i,4),Jf) max(ind(i,5),Ki) min(ind(i,6),Kf)];
    end
    toRemove = ind(:,3) > ind(:,4) | ind(:,5) > ind(:,6);
    ind(toRemove,:) = [];
    if size(vd,1) == 1 && n > 1 % If vd is a vector, this is a list of boundary conditions. The n > 1 condition checks is this is not a one-line matrix for when a single corner is present
        vd(toRemove(:)) = [];
    else
        vd(toRemove,:) = []; % If vd is a matrix, this is a list of corners
    end
    n = n-sum(toRemove);
    
    if strcmp(type,'n')
        for i = 1:n
            switch vd(i)
                case 3
                    if ind(i,4) + neumannLength > Jf
                        error('There is a y+ Neumann condition at J = %d crossing a domain slice at J = %d. Consider changing p_row.',ind(i,4),Jf)
                    end
                case 4
                    if ind(i,3) - neumannLength < Ji
                        error('There is a y- Neumann condition at J = %d crossing a domain slice at J = %d. Consider changing p_row.',ind(i,3),Ji)
                    end
                case 5
                    if ind(i,6) + neumannLength > Kf
                        error('There is a z+ Neumann condition at K = %d crossing a domain slice at K = %d. Consider changing p_col.',ind(i,6),Kf)
                    end
                case 6
                    if ind(i,5) - neumannLength < Ki
                        error('There is a z- Neumann condition at K = %d crossing a domain slice at K = %d. Consider changing p_col.',ind(i,5),Ki)
                    end
            end
        end
    end
end
