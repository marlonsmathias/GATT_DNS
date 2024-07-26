%% Define flow region
% Find which nodes will actually contain a flow and which ones will be in or at a wall
flowRegion = true(mesh.nx,mesh.ny,mesh.nz);

% Add cavities to the flow region
if isfield(flowType,'cav')
    for i = 1:length(flowType.cav)
        x = flowType.cav{i}.x;
        y = flowType.cav{i}.y;
        z = flowType.cav{i}.z;
        
        flowRegion(mesh.X>x(1) & mesh.X<x(2), mesh.Y>y(1) & mesh.Y<y(2), mesh.Z>z(1) & mesh.Z<z(2)) = true;
    end
end

% Remove roughnesses from the flow
if isfield(flowType,'rug')
    for i = 1:length(flowType.rug)
        x = flowType.rug{i}.x;
        y = flowType.rug{i}.y;
        z = flowType.rug{i}.z;
        
        flowRegion(mesh.X>=x(1) & mesh.X<=x(2), mesh.Y>=y(1) & mesh.Y<=y(2), mesh.Z>=z(1) & mesh.Z<=z(2)) = false;
    end
end

% Add outer walls
flowRegion(:,[1 end],:) = false;

%% Get walls
findWallsForBoundaries

%% Add walls to boundary conditions
for i = 1:6
    switch i
        case 1
            wallPosition = wallFrontLimits;
            wallDir = 'xi';
        case 2
            wallPosition = wallBackLimits;
            wallDir = 'xf';
        case 3
            wallPosition = wallUpLimits;
            wallDir = 'yi';
        case 4
            wallPosition = wallDownLimits;
            wallDir = 'yf';
        case 5
            wallPosition = wallRightLimits;
            wallDir = 'zi';
        case 6
            wallPosition = wallLeftLimits;
            wallDir = 'zf';
    end
    
    for j = 1:size(wallPosition,1)
        var{end+1} = 'p';  %#ok<*SAGROW>
        type{end+1} = 'neu';
        dir{end+1} = wallDir;
        val(end+1) = 0;
        xi(end+1) = wallPosition(j,1);
        xf(end+1) = wallPosition(j,2);
        yi(end+1) = wallPosition(j,3);
        yf(end+1) = wallPosition(j,4);
        zi(end+1) = wallPosition(j,5);
        zf(end+1) = wallPosition(j,6);
    end
        
end

%% Add regions that are inside walls
for i = 1:size(insideWalls,1)

    var{end+1} = 'p';
    type{end+1} = 'dir';
    dir{end+1} = 'yi';
    val(end+1) = P0;
    xi(end+1) = insideWalls(i,1);
    xf(end+1) = insideWalls(i,2);
    yi(end+1) = insideWalls(i,3);
    yf(end+1) = insideWalls(i,4);
    zi(end+1) = insideWalls(i,5);
    zf(end+1) = insideWalls(i,6);

    var{end+1} = 'u';
    type{end+1} = 'dir';
    dir{end+1} = 'yi';
    val(end+1) = 0;
    xi(end+1) = insideWalls(i,1);
    xf(end+1) = insideWalls(i,2);
    yi(end+1) = insideWalls(i,3);
    yf(end+1) = insideWalls(i,4);
    zi(end+1) = insideWalls(i,5);
    zf(end+1) = insideWalls(i,6);
    
    var{end+1} = 'v';
    type{end+1} = 'dir';
    dir{end+1} = 'yi';
    val(end+1) = 0;
    xi(end+1) = insideWalls(i,1);
    xf(end+1) = insideWalls(i,2);
    yi(end+1) = insideWalls(i,3);
    yf(end+1) = insideWalls(i,4);
    zi(end+1) = insideWalls(i,5);
    zf(end+1) = insideWalls(i,6);
    
    var{end+1} = 'w';
    type{end+1} = 'dir';
    dir{end+1} = 'yi';
    val(end+1) = 0;
    xi(end+1) = insideWalls(i,1);
    xf(end+1) = insideWalls(i,2);
    yi(end+1) = insideWalls(i,3);
    yf(end+1) = insideWalls(i,4);
    zi(end+1) = insideWalls(i,5);
    zf(end+1) = insideWalls(i,6);
    
	if isfield(flowType,'tWallRelative')
		eWall = E0 * flowType.tWallRelative;
	else
		eWall = E0;
	end
	
    var{end+1} = 'e';
    type{end+1} = 'dir';
    dir{end+1} = 'yi';
    val(end+1) = eWall;
    xi(end+1) = insideWalls(i,1);
    xf(end+1) = insideWalls(i,2);
    yi(end+1) = insideWalls(i,3);
    yf(end+1) = insideWalls(i,4);
    zi(end+1) = insideWalls(i,5);
    zf(end+1) = insideWalls(i,6);
end

% Add moving walls
var{end+1} = 'u';
type{end+1} = 'dir';
dir{end+1} = 'yi';
val(end+1) = flowParameters.lowerWallVelocity;
xi(end+1) = 1;
xf(end+1) = mesh.nx;
yi(end+1) = 1;
yf(end+1) = 1;
zi(end+1) = 1;
zf(end+1) = mesh.nz;

var{end+1} = 'u';
type{end+1} = 'dir';
dir{end+1} = 'yf';
val(end+1) = flowParameters.upperWallVelocity;
xi(end+1) = 1;
xf(end+1) = mesh.nx;
yi(end+1) = mesh.ny;
yf(end+1) = mesh.ny;
zi(end+1) = 1;
zf(end+1) = mesh.nz;