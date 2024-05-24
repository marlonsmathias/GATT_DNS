function [boundary,mesh] = getBoundaryConditions(flowType,mesh,flowParameters,neumannOrder)

% This function runs the relevant routine for the flowType defined in the parameters and outputs the variables in a structure
% The disturbances are also added to the structure
% The stencil and coefficients for Neumann conditions is defined here

%% Define base value for gamma, P and E
gamma = flowParameters.gamma;
E0 = 1/((gamma^2-gamma)*flowParameters.Ma^2);
P0 = (gamma-1)*E0;

%% Create variables that will recieve boundary information
var = {}; % Variable to which the condition is applied
type = {}; % Type of boundary. Either 'dir' for a Dirichlet condition or 'neu' for a Neumann condition
dir = {}; % Direction of the derivative for the Neumann condition (xi, xf, yi, yf, zi or zf)
val = []; % Value the Variable assumes in a Dirichlet condition or value of the derivative in a Neumann condition
xi = []; % Starting j index for the condition
xf = []; % Ending j index for the variable
yi = []; % Starting i index for the condition
yf = []; % Ending i index for the condition
zi = []; % Starting k index for the condition
zf = []; % Ending k index for the condition

%% Check the type of boundaries and call the appropriate subroutine
if exist(['source/boundaries/' flowType.name '.m'],'file')
    eval(flowType.name)
else
    error('Boundary condition file not found, check source/boundaries folder')
end

%% Join all relevant variables in the output structure
boundary.var = var;
boundary.type = type;
boundary.dir = dir;
boundary.val = val;
boundary.xi = xi;
boundary.xf = xf;
boundary.yi = yi;
boundary.yf = yf;
boundary.zi = zi;
boundary.zf = zf;

boundary.flowRegion = flowRegion;
boundary.E0 = E0;

boundary.wall.up = wallUpLimits;
boundary.wall.down = wallDownLimits;
boundary.wall.front = wallFrontLimits;
boundary.wall.back = wallBackLimits;
boundary.wall.right = wallRightLimits;
boundary.wall.left = wallLeftLimits;

boundary.corners = corners;

boundary.gamma = gamma;

%% Get wall region
boundary.insideWall = ~flowRegion;
for i = 1:6
    switch i
        case 1; currentWall = wallUpLimits;
        case 2; currentWall = wallDownLimits;
        case 3; currentWall = wallFrontLimits;
        case 4; currentWall = wallBackLimits;
        case 5; currentWall = wallRightLimits;
        case 6; currentWall = wallLeftLimits;
    end
    for j = 1:size(currentWall,1)
    boundary.insideWall(currentWall(j,1):currentWall(j,2),currentWall(j,3):currentWall(j,4),currentWall(j,5):currentWall(j,6)) = false;
    end
end

%% Get coefficients for neumann boundary conditions

NC1 = {[-1, 1], [-3/2, 2, -1/2], [-11/6, 3, -3/2, 1/3], [-25/12, 4, -3, 4/3, -1/4], [-137/60, 5, -5, 10/3, -5/4, 1/5], [-49/20, 6, -15/2, 20/3, -15/4, 6/5, -1/6]};

NC2 = {[1, -2, 1], [2, -5, 4, -1], [35/12, -26/3, 19/2, -14/3, 11/12], [15/4, -77/6, 107/6, -13, 61/12, -5/6], [203/45, -87/5, 117/4, -254/9, 33/2, -27/5, 137/180], [469/90, -223/10, 879/20, -949/18, 41, -201/10, 1019/180, -7/10]};

boundary.neumannCoeffs = -NC1{neumannOrder(1)}(2:end)/NC1{neumannOrder(1)}(1);
boundary.neumann2Coeffs = -NC2{neumannOrder(2)}(2:end)/NC2{neumannOrder(2)}(1);

%% Add all disturbances to boundary structure
boundary.disturb = {};
if isfield(flowType,'disturb')
    for i = 1:length(flowType.disturb)
        if flowType.disturb{i}.active
            boundary.disturb{i}.type = flowType.disturb{i}.type;

            if isfield(flowType.disturb{i},'forcing') && flowType.disturb{i}.forcing
                boundary.disturb{i}.forcing=true;
            else
                boundary.disturb{i}.forcing=false;
            end
            
            if isfield(flowType.disturb{i},'par')
				if isnumeric(flowType.disturb{i}.par) % Convert vector to cell, if needed
					flowType.disturb{i}.par = num2cell(flowType.disturb{i}.par);
				end
                boundary.disturb{i}.par = flowType.disturb{i}.par;
            else
                boundary.disturb{i}.par = {};
            end
            
            boundary.disturb{i}.var = flowType.disturb{i}.var;
            
            xi = flowType.disturb{i}.x(1);
            xf = flowType.disturb{i}.x(2);
			
			xi(isinf(xi)&&xi<0) = mesh.X(1);
            xf(isinf(xf)&&xf<0) = mesh.X(1);
			xi(isinf(xi)&&xi>0) = mesh.X(end);
            xf(isinf(xf)&&xf>0) = mesh.X(end);
			
            yi = flowType.disturb{i}.y(1);
            yf = flowType.disturb{i}.y(2);

			yi(isinf(yi)&&yi<0) = mesh.Y(1);
            yf(isinf(yf)&&yf<0) = mesh.Y(1);
			yi(isinf(yi)&&yi>0) = mesh.Y(end);
            yf(isinf(yf)&&yf>0) = mesh.Y(end);
			
            zi = flowType.disturb{i}.z(1);
            zf = flowType.disturb{i}.z(2);
            
			zi(isinf(zi)&&zi<0) = mesh.Z(1);
            zf(isinf(zf)&&zf<0) = mesh.Z(1);
			zi(isinf(zi)&&zi>0) = mesh.Z(end);
            zf(isinf(zf)&&zf>0) = mesh.Z(end);
            
            if ~isfield(flowType.disturb{i},'extraNodes')
                boundary.disturb{i}.extraNodes = [0 0 0 0 0 0];
            else
                boundary.disturb{i}.extraNodes = flowType.disturb{i}.extraNodes;
            end
            
            [~,boundary.disturb{i}.ind(1)] = find(mesh.X>=xi,1,'first');
            [~,boundary.disturb{i}.ind(2)] = find(mesh.X<=xf,1,'last');
            [~,boundary.disturb{i}.ind(3)] = find(mesh.Y>=yi,1,'first');
            [~,boundary.disturb{i}.ind(4)] = find(mesh.Y<=yf,1,'last');
            [~,boundary.disturb{i}.ind(5)] = find(mesh.Z>=zi,1,'first');
            [~,boundary.disturb{i}.ind(6)] = find(mesh.Z<=zf,1,'last');
            boundary.disturb{i}.ind(1) = boundary.disturb{i}.ind(1) - boundary.disturb{i}.extraNodes(1);
            boundary.disturb{i}.ind(2) = boundary.disturb{i}.ind(2) + boundary.disturb{i}.extraNodes(2);
            boundary.disturb{i}.ind(3) = boundary.disturb{i}.ind(3) - boundary.disturb{i}.extraNodes(3);
            boundary.disturb{i}.ind(4) = boundary.disturb{i}.ind(4) + boundary.disturb{i}.extraNodes(4);
            boundary.disturb{i}.ind(5) = boundary.disturb{i}.ind(5) - boundary.disturb{i}.extraNodes(5);
            boundary.disturb{i}.ind(6) = boundary.disturb{i}.ind(6) + boundary.disturb{i}.extraNodes(6);
        end
    end
end

end
