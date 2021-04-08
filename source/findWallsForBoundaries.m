% This script is called by getBoundaryConditions.m, it takes flowRegion as input and outputs a list of walls and corners in the domain
% Each wall is defined by six columns in a list. The columns stand for [xi xf yi yf zi zf]
% Corners are defined by their indices (corners.limits) and they directions (corners.dir). For example, a corner defined by direction [-1 1 0] faces backwards in x, upwards in y and is constant in z
% corners.adiabatic indicates if the corner is adiabatic or isothermal. It defaults to isothermal.
% Regions that are completely confined inside walls are placed in the insideWalls variable

%% Find and remove infinitelly thin walls
% in x
i = 2:mesh.nx-1;
flowRegion(i,:,:) =  flowRegion(i,:,:) | (flowRegion(i-1,:,:) & flowRegion(i+1,:,:));

% in y
j = 2:mesh.ny-1;
flowRegion(:,j,:) =  flowRegion(:,j,:) | (flowRegion(:,j-1,:) & flowRegion(:,j+1,:));

% in z
k = 2:mesh.nz-1;
flowRegion(:,:,k) =  flowRegion(:,:,k) | (flowRegion(:,:,k-1) & flowRegion(:,:,k+1));

%% Find walls
wallFront = false(mesh.nx,mesh.ny,mesh.nz);
wallFront(1:end-1,:,:) = diff(flowRegion,[],1) == 1;
wallBack = false(mesh.nx,mesh.ny,mesh.nz);
wallBack(2:end,:,:) = diff(flowRegion,[],1) == -1;

wallUp = false(mesh.nx,mesh.ny,mesh.nz);
wallUp(:,1:end-1,:) = diff(flowRegion,[],2) == 1;
wallDown = false(mesh.nx,mesh.ny,mesh.nz);
wallDown(:,2:end,:) = diff(flowRegion,[],2) == -1;

wallRight = false(mesh.nx,mesh.ny,mesh.nz);
wallRight(:,:,1:end-1) = diff(flowRegion,[],3) == 1;
wallLeft = false(mesh.nx,mesh.ny,mesh.nz);
wallLeft(:,:,2:end) = diff(flowRegion,[],3) == -1;

%% Find regions for boundaries in walls
% In x
wallFrontLimits = [];
wallBackLimits = [];

for i = 1:mesh.nx
    localWall = wallFront(i,:,:); % Find front facing walls
    if any(localWall(:))
        for k = 1:mesh.nz
            if k == 1 || isempty(wallStarts) % Get starting and ending points in y
                wallStarts = find([localWall(1,1,k) squeeze(diff(localWall(1,:,k),[],2)==1)]);
                wallEnds = find([squeeze(diff(localWall(1,:,k),[],2)==-1) localWall(1,end,k)]);
                kStart = k;
            end
            if k == mesh.nz || any(localWall(1,:,k) ~= localWall(1,:,k+1)) % If reached end in z or wall change, record boundary
                for j = 1:length(wallStarts)
                    wallFrontLimits(end+1,:) = [i i wallStarts(j) wallEnds(j) kStart k]; %#ok<*SAGROW>
                end
                wallStarts = [];
                wallEnds = [];
            end
        end
    end
    
    localWall = wallBack(i,:,:); % Find down facing walls
    if any(localWall(:))
        for k = 1:mesh.nz
            if k == 1 || isempty(wallStarts) % Get starting and ending points in y
                wallStarts = find([localWall(1,1,k) squeeze(diff(localWall(1,:,k),[],2)==1)]);
                wallEnds = find([squeeze(diff(localWall(1,:,k),[],2)==-1) localWall(1,end,k)]);
                kStart = k;
            end
            if k == mesh.nz || any(localWall(1,:,k) ~= localWall(1,:,k+1)) % If reached end in z or wall change, record boundary
                for j = 1:length(wallStarts)
                    wallBackLimits(end+1,:) = [i i wallStarts(j) wallEnds(j) kStart k];
                end
                wallStarts = [];
                wallEnds = [];
            end
        end
    end
end

% in y
wallUpLimits = [];
wallDownLimits = [];

for j = 1:mesh.ny
    localWall = wallUp(:,j,:); % Find front facing walls
    if any(localWall(:))
        for k = 1:mesh.nz
            if k == 1 || isempty(wallStarts) % Get starting and ending points in x
                wallStarts = find([localWall(1,1,k) squeeze(diff(localWall(:,1,k),[],1)==1)']);
                wallEnds = find([squeeze(diff(localWall(:,1,k),[],1)==-1)' localWall(end,1,k)]);
                kStart = k;
            end
            if k == mesh.nz || any(localWall(:,1,k) ~= localWall(:,1,k+1)) % If reached end in z or wall change, record boundary
                for i = 1:length(wallStarts)
                    wallUpLimits(end+1,:) = [wallStarts(i) wallEnds(i) j j kStart k];
                end
                wallStarts = [];
                wallEnds = [];
            end
        end
    end
    
    localWall = wallDown(:,j,:); % Find front facing walls
    if any(localWall(:))
        for k = 1:mesh.nz
            if k == 1 || isempty(wallStarts) % Get starting and ending points in x
                wallStarts = find([localWall(1,1,k) squeeze(diff(localWall(:,1,k),[],1)==1)']);
                wallEnds = find([squeeze(diff(localWall(:,1,k),[],1)==-1)' localWall(end,1,k)]);
                kStart = k;
            end
            if k == mesh.nz || any(localWall(:,1,k) ~= localWall(:,1,k+1)) % If reached end in z or wall change, record boundary
                for i = 1:length(wallStarts)
                    wallDownLimits(end+1,:) = [wallStarts(i) wallEnds(i) j j kStart k];
                end
                wallStarts = [];
                wallEnds = [];
            end
        end
    end
end

% in z
wallRightLimits = [];
wallLeftLimits = [];

for k = 1:mesh.nz
    localWall = wallRight(:,:,k); % Find left facing walls
    if any(localWall(:))
        for i = 1:mesh.nx
            if i == 1 || isempty(wallStarts) % Get starting and ending points in z
                wallStarts = find([localWall(i,1,1) squeeze(diff(localWall(i,:,1),[],2)==1)]);
                wallEnds = find([squeeze(diff(localWall(i,:,1),[],2)==-1) localWall(i,end,1)]);
                iStart = i;
            end
            if i == mesh.nx || any(localWall(i,:,1) ~= localWall(i+1,:,1)) % If reached end in x or wall change, record boundary
                for j = 1:length(wallStarts)
                    wallRightLimits(end+1,:) = [iStart i wallStarts(j) wallEnds(j) k k];
                end
                wallStarts = [];
                wallEnds = [];
            end
        end
    end

    localWall = wallLeft(:,:,k); % Find right facing walls
    if any(localWall(:))
        for i = 1:mesh.nx
            if i == 1 || isempty(wallStarts) % Get starting and ending points in z
                wallStarts = find([localWall(i,1,1) squeeze(diff(localWall(i,:,1),[],2)==1)]);
                wallEnds = find([squeeze(diff(localWall(i,:,1),[],2)==-1) localWall(i,end,1)]);
                iStart = i;
            end
            if i == mesh.nx || any(localWall(i,:,1) ~= localWall(i+1,:,1)) % If reached end in x or wall change, record boundary
                for j = 1:length(wallStarts)
                    wallLeftLimits(end+1,:) = [iStart i wallStarts(j) wallEnds(j) k k];
                end
                wallStarts = [];
                wallEnds = [];
            end
        end
    end
    
end

%% Join adjacent walls

vars = {'wallFrontLimits','wallBackLimits','wallUpLimits','wallDownLimits'};
for n = 1:4
    done = false;
    eval(['currentWall = ' vars{n} ';']);
    while ~done
        done = true;
        nWalls = size(currentWall,1);
        for i = 1:nWalls-1
            for j = i+1:nWalls
                if all(currentWall(i,1:4) == currentWall(j,1:4)) && currentWall(i,6) == currentWall(j,5)-1
                    toMerge = [i j];
                    done = false;
                end
            end
            if done == false
                break
            end
        end
        if done == false
            currentWall(toMerge(1),6) = currentWall(toMerge(2),6);
            currentWall(toMerge(2),:) = [];
        end
    end
    eval([vars{n} ' = currentWall;']);
end

%% Find corners

corners.limits = [];
corners.dir = [];

% With 2 walls

% Constant z
cornersMatrix(:,:,:,1) = wallFront & wallUp;
cornersMatrix(:,:,:,2) = wallFront & wallDown;
cornersMatrix(:,:,:,3) = wallBack & wallUp;
cornersMatrix(:,:,:,4) = wallBack & wallDown;
cornerDirections = [1 1 0; 1 -1 0; -1 1 0; -1 -1 0];

for i = 1:mesh.nx
    for j = 1:mesh.ny
        for m = 1:4
            cornerRow = squeeze(cornersMatrix(i,j,:,m))';
            if any(cornerRow)
                cornerStarts = find([cornerRow(1) diff(cornerRow)==1]);
                cornerEnds = find([diff(cornerRow)==-1 cornerRow(end)]);
                for n = 1:length(cornerStarts)
                    corners.limits(end+1,:) = [i i j j cornerStarts(n) cornerEnds(n)];
                    corners.dir(end+1,:) = cornerDirections(m,:);
                end
            end
        end
    end
end

% Constant x
cornersMatrix(:,:,:,1) = wallFront & wallRight;
cornersMatrix(:,:,:,2) = wallFront & wallLeft;
cornersMatrix(:,:,:,3) = wallBack & wallRight;
cornersMatrix(:,:,:,4) = wallBack & wallLeft;
cornerDirections = [1 0 1; 1 0 -1; -1 0 1; -1 0 -1];

for i = 1:mesh.nx
    for k = 1:mesh.nz
        for m = 1:4
            cornerRow = cornersMatrix(i,:,k,m);
            if any(cornerRow)
                cornerStarts = find([cornerRow(1) diff(cornerRow)==1]);
                cornerEnds = find([diff(cornerRow)==-1 cornerRow(end)]);
                for n = 1:length(cornerStarts)
                    corners.limits(end+1,:) = [i i cornerStarts(n) cornerEnds(n) k k];
                    corners.dir(end+1,:) = cornerDirections(m,:);
                end
            end
        end
    end
end

% Constant y
cornersMatrix(:,:,:,1) = wallUp & wallRight;
cornersMatrix(:,:,:,2) = wallUp & wallLeft;
cornersMatrix(:,:,:,3) = wallDown & wallRight;
cornersMatrix(:,:,:,4) = wallDown & wallLeft;
cornerDirections = [0 1 1; 0 1 -1; 0 -1 1; 0 -1 -1];

for j = 1:mesh.ny
    for k = 1:mesh.nz
        for m = 1:4
            cornerRow = cornersMatrix(:,j,k,m)';
            if any(cornerRow)
                cornerStarts = find([cornerRow(1) diff(cornerRow)==1]);
                cornerEnds = find([diff(cornerRow)==-1 cornerRow(end)]);
                for n = 1:length(cornerStarts)
                    corners.limits(end+1,:) = [cornerStarts(n) cornerEnds(n) j j k k];
                    corners.dir(end+1,:) = cornerDirections(m,:);
                end
            end
        end
    end
end

% With 3 walls
cornersMatrix(:,:,:,1) = wallFront & wallUp & wallRight;
cornersMatrix(:,:,:,2) = wallFront & wallUp & wallLeft;
cornersMatrix(:,:,:,3) = wallFront & wallDown & wallRight;
cornersMatrix(:,:,:,4) = wallFront & wallDown & wallLeft;
cornersMatrix(:,:,:,5) = wallBack & wallUp & wallRight;
cornersMatrix(:,:,:,6) = wallBack & wallUp & wallLeft;
cornersMatrix(:,:,:,7) = wallBack & wallDown & wallRight;
cornersMatrix(:,:,:,8) = wallBack & wallDown & wallLeft;

cornerDirections = [1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1; -1 1 1; -1 1 -1; -1 -1 1; -1 -1 -1];

[I,J,K,M] = ind2sub([mesh.nx,mesh.ny,mesh.nz,8],find(cornersMatrix(:,:,:,:)));
corners.limits = [corners.limits; I I J J K K];
corners.dir = [corners.dir; cornerDirections(M,:)];

corners.adiabatic = zeros(size(corners.dir,1),1);

%% Find regions completely contained inside walls
isWall = ~flowRegion;

% Find wall limits in x

insideWalls = [];
for k = 1:mesh.nz
    for j = 1:mesh.ny
        wallBegins = find([isWall(1,j,k); diff(isWall(:,j,k))==1]);
        wallEnds = find([diff(isWall(:,j,k))==-1; isWall(end,j,k)]);
        for i = 1:length(wallBegins)
            insideWalls(end+1,:) = [wallBegins(i) wallEnds(i) j j k k];
        end
    end
end

% Sort by same wall limits
if ~isempty(insideWalls)
    [~,~,ind] = unique(insideWalls(:,[1 2]),'rows','stable');
    [~,ind] = sort(ind);
    insideWalls = insideWalls(ind,:);
end

% Merge limits in y
i = 1;
while i < size(insideWalls,1)
    if all(insideWalls(i,[1 2 5 6]) == insideWalls(i+1,[1 2 5 6])) && insideWalls(i,4)+1 == insideWalls(i+1,3)
        insideWalls(i,4) = insideWalls(i+1,4);
        insideWalls(i+1,:) = [];
    else
        i = i + 1;
    end
end

% Merge limits in z
i = 1;
while i < size(insideWalls,1)
    if all(insideWalls(i,[1 2 3 4]) == insideWalls(i+1,[1 2 3 4])) && insideWalls(i,6)+1 == insideWalls(i+1,5)
        insideWalls(i,6) = insideWalls(i+1,6);
        insideWalls(i+1,:) = [];
    else
        i = i + 1;
    end
end
