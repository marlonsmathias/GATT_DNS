%% Subroutine for an isothermal boundary layer with or without cavities and roughtnesses
% To be called from getBoundaryConditions

%% Inflow
if mesh.X(1) <= 0
	% u
	var{end+1} = 'u';
	type{end+1} = 'dir';
	dir{end+1} = 'xi';
	val(end+1) = 1;
	xi(end+1) = 1;
	xf(end+1) = 1;
	yi(end+1) = 1;
	yf(end+1) = mesh.ny;
	zi(end+1) = 1;
	zf(end+1) = mesh.nz;

	% v
	var{end+1} = 'v';
	type{end+1} = 'dir';
	dir{end+1} = 'xi';
	val(end+1) = 0;
	xi(end+1) = 1;
	xf(end+1) = 1;
	yi(end+1) = 1;
	yf(end+1) = mesh.ny;
	zi(end+1) = 1;
	zf(end+1) = mesh.nz;

	% w
	var{end+1} = 'w';
	type{end+1} = 'dir';
	dir{end+1} = 'xi';
	val(end+1) = 0;
	xi(end+1) = 1;
	xf(end+1) = 1;
	yi(end+1) = 1;
	yf(end+1) = mesh.ny;
	zi(end+1) = 1;
	zf(end+1) = mesh.nz;

	% e
	var{end+1} = 'e';
	type{end+1} = 'dir';
	dir{end+1} = 'xi';
	val(end+1) = E0;
	xi(end+1) = 1;
	xf(end+1) = 1;
	yi(end+1) = 1;
	yf(end+1) = mesh.ny;
	zi(end+1) = 1;
	zf(end+1) = mesh.nz;
else
	if ~isfield(flowType,'disturb')
		flowType.disturb = cell(0);
	else
		flowType.disturb(2:end+1) = flowType.disturb;
	end
	flowType.disturb{1} = [];
	flowType.disturb{1}.x = [mesh.X(1) mesh.X(1)];
	flowType.disturb{1}.y = [-inf inf];
	flowType.disturb{1}.z = [-inf inf];
	flowType.disturb{1}.var = 'UVRWE';
	flowType.disturb{1}.type = 'holdInlet';
	flowType.disturb{1}.active = true;
	flowType.disturb{1}.par = [0]; %Hold density
end

% p
var{end+1} = 'p';
type{end+1} = 'neu';
dir{end+1} = 'xi';
val(end+1) = 0;
xi(end+1) = 1;
xf(end+1) = 1;
yi(end+1) = 1;
yf(end+1) = mesh.ny;
zi(end+1) = 1;
zf(end+1) = mesh.nz;


%% Outflow
% u
var{end+1} = 'u';
type{end+1} = 'sec';
dir{end+1} = 'xf';
val(end+1) = 0;
xi(end+1) = mesh.nx;
xf(end+1) = mesh.nx;
yi(end+1) = 1;
yf(end+1) = mesh.ny;
zi(end+1) = 1;
zf(end+1) = mesh.nz;

% v
var{end+1} = 'v';
type{end+1} = 'sec';
dir{end+1} = 'xf';
val(end+1) = 0;
xi(end+1) = mesh.nx;
xf(end+1) = mesh.nx;
yi(end+1) = 1;
yf(end+1) = mesh.ny;
zi(end+1) = 1;
zf(end+1) = mesh.nz;

% w
var{end+1} = 'w';
type{end+1} = 'sec';
dir{end+1} = 'xf';
val(end+1) = 0;
xi(end+1) = mesh.nx;
xf(end+1) = mesh.nx;
yi(end+1) = 1;
yf(end+1) = mesh.ny;
zi(end+1) = 1;
zf(end+1) = mesh.nz;

% e
var{end+1} = 'e';
type{end+1} = 'sec';
dir{end+1} = 'xf';
val(end+1) = 0;
xi(end+1) = mesh.nx;
xf(end+1) = mesh.nx;
yi(end+1) = 1;
yf(end+1) = mesh.ny;
zi(end+1) = 1;
zf(end+1) = mesh.nz;

% p
var{end+1} = 'p';
type{end+1} = 'dir';
dir{end+1} = 'xf';
val(end+1) = P0;
xi(end+1) = mesh.nx;
xf(end+1) = mesh.nx;
yi(end+1) = 1;
yf(end+1) = mesh.ny;
zi(end+1) = 1;
zf(end+1) = mesh.nz;


%% Outerflow
% u
var{end+1} = 'u';
type{end+1} = 'sec';
dir{end+1} = 'yf';
val(end+1) = 0;
xi(end+1) = 1;
xf(end+1) = mesh.nx;
yi(end+1) = mesh.ny;
yf(end+1) = mesh.ny;
zi(end+1) = 1;
zf(end+1) = mesh.nz;

% v
var{end+1} = 'v';
type{end+1} = 'sec';
dir{end+1} = 'yf';
val(end+1) = 0;
xi(end+1) = 1;
xf(end+1) = mesh.nx;
yi(end+1) = mesh.ny;
yf(end+1) = mesh.ny;
zi(end+1) = 1;
zf(end+1) = mesh.nz;

% w
var{end+1} = 'w';
type{end+1} = 'sec';
dir{end+1} = 'yf';
val(end+1) = 0;
xi(end+1) = 1;
xf(end+1) = mesh.nx;
yi(end+1) = mesh.ny;
yf(end+1) = mesh.ny;
zi(end+1) = 1;
zf(end+1) = mesh.nz;

% e
var{end+1} = 'e';
type{end+1} = 'sec';
dir{end+1} = 'yf';
val(end+1) = 0;
xi(end+1) = 1;
xf(end+1) = mesh.nx;
yi(end+1) = mesh.ny;
yf(end+1) = mesh.ny;
zi(end+1) = 1;
zf(end+1) = mesh.nz;

% p
var{end+1} = 'p';
type{end+1} = 'sec';
dir{end+1} = 'yf';
val(end+1) = 0;
xi(end+1) = 1;
xf(end+1) = mesh.nx;
yi(end+1) = mesh.ny;
yf(end+1) = mesh.ny;
zi(end+1) = 1;
zf(end+1) = mesh.nz;

%% Outerflow for non-periodic 3D
if mesh.nz > 1 && ~mesh.z.periodic

    % u
    var{end+1} = 'u';
    type{end+1} = 'neu';
    dir{end+1} = 'zi';
    val(end+1) = 0;
    xi(end+1) = 1;
    xf(end+1) = mesh.nx;
    yi(end+1) = 1;
    yf(end+1) = mesh.ny;
    zi(end+1) = 1;
    zf(end+1) = 1;

    % v
    var{end+1} = 'v';
    type{end+1} = 'neu';
    dir{end+1} = 'zi';
    val(end+1) = 0;
    xi(end+1) = 1;
    xf(end+1) = mesh.nx;
    yi(end+1) = 1;
    yf(end+1) = mesh.ny;
    zi(end+1) = 1;
    zf(end+1) = 1;

    % w
    var{end+1} = 'w';
    type{end+1} = 'neu';
    dir{end+1} = 'zi';
    val(end+1) = 0;
    xi(end+1) = 1;
    xf(end+1) = mesh.nx;
    yi(end+1) = 1;
    yf(end+1) = mesh.ny;
    zi(end+1) = 1;
    zf(end+1) = 1;

    % e
    var{end+1} = 'e';
    type{end+1} = 'neu';
    dir{end+1} = 'zi';
    val(end+1) = 0;
    xi(end+1) = 1;
    xf(end+1) = mesh.nx;
    yi(end+1) = 1;
    yf(end+1) = mesh.ny;
    zi(end+1) = 1;
    zf(end+1) = 1;

    % p
    var{end+1} = 'p';
    type{end+1} = 'neu';
    dir{end+1} = 'zi';
    val(end+1) = 0;
    xi(end+1) = 1;
    xf(end+1) = mesh.nx;
    yi(end+1) = 1;
    yf(end+1) = mesh.ny;
    zi(end+1) = 1;
    zf(end+1) = 1;
    
    % u
    var{end+1} = 'u';
    type{end+1} = 'neu';
    dir{end+1} = 'zf';
    val(end+1) = 0;
    xi(end+1) = 1;
    xf(end+1) = mesh.nx;
    yi(end+1) = 1;
    yf(end+1) = mesh.ny;
    zi(end+1) = mesh.nz;
    zf(end+1) = mesh.nz;

    % v
    var{end+1} = 'v';
    type{end+1} = 'neu';
    dir{end+1} = 'zf';
    val(end+1) = 0;
    xi(end+1) = 1;
    xf(end+1) = mesh.nx;
    yi(end+1) = 1;
    yf(end+1) = mesh.ny;
    zi(end+1) = mesh.nz;
    zf(end+1) = mesh.nz;

    % w
    var{end+1} = 'w';
    type{end+1} = 'neu';
    dir{end+1} = 'zf';
    val(end+1) = 0;
    xi(end+1) = 1;
    xf(end+1) = mesh.nx;
    yi(end+1) = 1;
    yf(end+1) = mesh.ny;
    zi(end+1) = mesh.nz;
    zf(end+1) = mesh.nz;

    % e
    var{end+1} = 'e';
    type{end+1} = 'neu';
    dir{end+1} = 'zf';
    val(end+1) = 0;
    xi(end+1) = 1;
    xf(end+1) = mesh.nx;
    yi(end+1) = 1;
    yf(end+1) = mesh.ny;
    zi(end+1) = mesh.nz;
    zf(end+1) = mesh.nz;

    % p
    var{end+1} = 'p';
    type{end+1} = 'neu';
    dir{end+1} = 'zf';
    val(end+1) = 0;
    xi(end+1) = 1;
    xf(end+1) = mesh.nx;
    yi(end+1) = 1;
    yf(end+1) = mesh.ny;
    zi(end+1) = mesh.nz;
    zf(end+1) = mesh.nz;

end

%% Define flow region
% Find which nodes will actually contain a flow and which ones will be in or at a wall
flowRegion = true(mesh.nx,mesh.ny,mesh.nz);

% Add flat plate
[~,wallJ] = min(abs(mesh.Y));
flowRegion(:,1:wallJ,:) = false;

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

%% Get walls
findWallsForBoundaries

%% Create wall for free-slip region
if mesh.X(1) < 0
    wallUpLimits(:,1) = max(wallUpLimits(:,1),find(mesh.X>=0,1,'first')); % Move existing walls
    wallUpLimits(end+1,:) = [1 find(mesh.X<0,1,'last') [1 1]*find(mesh.Y>=0,1,'first') 1 mesh.nz];
    
    mesh.x.breakPoint = [find(mesh.X<0,1,'last') [1 1]*find(mesh.Y>=0,1,'first') 1 mesh.nz];
end

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

%% Add free slip wall, it is the last up facing wall
if mesh.X(1) < 0
    var{end+1} = 'u';
    type{end+1} = 'neu';
    dir{end+1} = 'yi';
    val(end+1) = 0;
    xi(end+1) = wallUpLimits(end,1);
    xf(end+1) = wallUpLimits(end,2);
    yi(end+1) = wallUpLimits(end,3);
    yf(end+1) = wallUpLimits(end,4);
    zi(end+1) = wallUpLimits(end,5);
    zf(end+1) = wallUpLimits(end,6);

    var{end+1} = 'w';
    type{end+1} = 'neu';
    dir{end+1} = 'yi';
    val(end+1) = 0;
    xi(end+1) = wallUpLimits(end,1);
    xf(end+1) = wallUpLimits(end,2);
    yi(end+1) = wallUpLimits(end,3);
    yf(end+1) = wallUpLimits(end,4);
    zi(end+1) = wallUpLimits(end,5);
    zf(end+1) = wallUpLimits(end,6);
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
