% This script is called by makeMatrices.m and gets the flowRegion array as input and outputs three arrays, one for each direction
% For example, typeMapX is ny by nz and each entry contains the type of the derivative to be used in the row
% derivStartsX and derivEndsX contain the begining and ending of each derivative type in typeMapX

%% In X
% Get matrix
flowRegion = permute(boundary.flowRegion,[2 1 3]);

% Allocate variables
derivStartsX = {};
derivEndsX = {};

% Find all different types of geometries
C = [];
for k = 1:mesh.nz
    C = [C; unique(flowRegion(:,:,k),'rows','sorted')]; %#ok<*AGROW>
end
C = unique(C,'rows','sorted');

% Identify where these regions are present
typeMapX = zeros(mesh.ny,mesh.nz);
for j = 1:mesh.ny
    for k = 1:mesh.nz
        typeMapX(j,k) = find(all(bsxfun(@eq,flowRegion(j,:,k),C),2));
    end
end

for i = 1:size(C,1) % Find walls in each of the regions types
    if any(C(i,:) ~= 0) % Only record regions with flow
        if mesh.x.periodic
            derivStartsX{end+1} = unique(find(diff(C(i,:))==1)); %#ok<*SAGROW>
            derivEndsX{end+1} = unique(find(diff(C(i,:))==-1)+1);
        else
            derivStartsX{end+1} = unique([1 find(diff(C(i,:))==1)]);
            derivEndsX{end+1} = unique([find(diff(C(i,:))==-1)+1 mesh.nx]);
        end
    else
        typeMapX(typeMapX == i) = 0;
        typeMapX(typeMapX > i) = typeMapX(typeMapX > 1) - 1;
    end
end

typeMapX(typeMapX == 0) = max(typeMapX(:));

if isfield(mesh.x,'breakPoint')
    for j = 1:size(mesh.x.breakPoint,1)
        ind = mesh.x.breakPoint(j,2:5);
        types = unique(typeMapX(ind(1):ind(2),ind(3):ind(4)));
        for i = 1:length(types)
            typeN = max(typeMapX(:))+1;
            derivStartsX{typeN} = derivStartsX{types(i)};
            derivEndsX{typeN} = derivEndsX{types(i)};

            derivStartsX{typeN}(end+1) = mesh.x.breakPoint(j,1)+1;
            derivEndsX{typeN}(end+1) = mesh.x.breakPoint(j,1);
            
            typeMapTemp = typeMapX(ind(1):ind(2),ind(3):ind(4));
            typeMapTemp(typeMapTemp == types(i)) = typeN;
            typeMapX(ind(1):ind(2),ind(3):ind(4)) = typeMapTemp;
        end
    end
end

%% In Y
% Get matrix
flowRegion = boundary.flowRegion;

% Allocate variables
derivStartsY = {};
derivEndsY = {};

% Find all different types of geometries
C = [];
for k = 1:mesh.nz
    C = [C; unique(flowRegion(:,:,k),'rows','sorted')];
end
C = unique(C,'rows','sorted');

% Identify where these regions are present
typeMapY = zeros(mesh.nx,mesh.nz);
for i = 1:mesh.nx
    for k = 1:mesh.nz
        typeMapY(i,k) = find(all(bsxfun(@eq,flowRegion(i,:,k),C),2));
    end
end

for i = 1:size(C,1) % Find walls in each of the regions types
    if any(C(i,:) ~= 0) % Only record regions with flow
        if mesh.y.periodic
            derivStartsY{end+1} = unique(find(diff(C(i,:))==1));
            derivEndsY{end+1} = unique(find(diff(C(i,:))==-1)+1);
        else
            derivStartsY{end+1} = unique([1 find(diff(C(i,:))==1)]);
            derivEndsY{end+1} = unique([find(diff(C(i,:))==-1)+1 mesh.ny]);
        end
    else
        typeMapY(typeMapY == i) = 0;
        typeMapY(typeMapY > i) = typeMapY(typeMapY > 1) - 1;
    end
end

typeMapY(typeMapY == 0) = max(typeMapY(:));

if isfield(mesh.y,'breakPoint')
    for j = 1:size(mesh.y.breakPoint,1)
        ind = mesh.y.breakPoint(j,2:5);
        types = unique(typeMapY(ind(1):ind(2),ind(3):ind(4)));
        for i = 1:length(types)
            typeN = max(typeMapY(:))+1;
            derivStartsY{typeN} = derivStartsY{types(i)};
            derivEndsY{typeN} = derivEndsY{types(i)};

            derivStartsY{typeN}(end+1) = mesh.y.breakPoint(j,1)+1;
            derivEndsY{typeN}(end+1) = mesh.y.breakPoint(j,1);
            
            typeMapTemp = typeMapY(ind(1):ind(2),ind(3):ind(4));
            typeMapTemp(typeMapTemp == types(i)) = typeN;
            typeMapY(ind(1):ind(2),ind(3):ind(4)) = typeMapTemp;
        end
    end
end

%% In Z
% Get matrix
flowRegion = permute(boundary.flowRegion,[1 3 2]);

% Allocate variables
derivStartsZ = {};
derivEndsZ = {};

% Find all different types of geometries
C = [];
for j = 1:mesh.ny
    C = [C; unique(flowRegion(:,:,j),'rows','sorted')];
end
C = unique(C,'rows','sorted');

% Identify where these regions are present
typeMapZ = zeros(mesh.nx,mesh.ny);
for i = 1:mesh.nx
    for j = 1:mesh.ny
        typeMapZ(i,j) = find(all(bsxfun(@eq,flowRegion(i,:,j),C),2));
    end
end

for i = 1:size(C,1) % Find walls in each of the regions types
    if any(C(i,:) ~= 0) % Only record regions with flow
        if mesh.z.periodic
            derivStartsZ{end+1} = unique(find(diff(C(i,:))==1));
            derivEndsZ{end+1} = unique(find(diff(C(i,:))==-1)+1);
        else
            derivStartsZ{end+1} = unique([1 find(diff(C(i,:))==1)]);
            derivEndsZ{end+1} = unique([find(diff(C(i,:))==-1)+1 mesh.nz]);
        end
    else
        typeMapZ(typeMapZ == i) = 0;
        typeMapZ(typeMapZ > i) = typeMapZ(typeMapZ > 1) - 1;
    end
end

typeMapZ(typeMapZ == 0) = max(typeMapZ(:));

if isfield(mesh.z,'breakPoint')
    for j = 1:size(mesh.z.breakPoint,1)
        ind = mesh.z.breakPoint(j,2:5);
        types = unique(typeMapZ(ind(1):ind(2),ind(3):ind(4)));
        for i = 1:length(types)
            typeN = max(typeMapZ(:))+1;
            derivStartsZ{typeN} = derivStartsZ{types(i)};
            derivEndsZ{typeN} = derivEndsZ{types(i)};

            derivStartsZ{typeN}(end+1) = mesh.z.breakPoint(j,1)+1;
            derivEndsZ{typeN}(end+1) = mesh.z.breakPoint(j,1);
            
            typeMapTemp = typeMapZ(ind(1):ind(2),ind(3):ind(4));
            typeMapTemp(typeMapTemp == types(i)) = typeN;
            typeMapZ(ind(1):ind(2),ind(3):ind(4)) = typeMapTemp;
        end
    end
end