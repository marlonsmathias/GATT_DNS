delta = 1;
Re = 1000;
Ma = 0.1;
Lx = 1;
Ly = Lx;
Lz = Lx;
nx = 64;
ny = nx;
nz = nx;
epsilon = 1;

%% Case name
caseName = ['TurbulenceBox-Re' num2str(Re) '-Ma' strrep(num2str(Ma),'.','') '-eps' num2str(epsilon) '-L' num2str(Lx) 'x' num2str(Ly) 'x' num2str(Lz) '-n' num2str(nx) 'x' num2str(ny) 'x' num2str(nz)];


%% Domain decomposition
p_row = 4;
p_col = 4;

%% Flow parameters
flowParameters.Re = Re;
flowParameters.Ma = Ma;
flowParameters.Pr = 0.71;
flowParameters.gamma = 1.4;
flowParameters.T0 = 300;

%% Domain parameters
domain.xi = 0;
domain.xf = Lx;
domain.yi = 0;
domain.yf = Ly;
domain.zi = 0;
domain.zf = Lz;

%% Flow type
flowType.name = 'periodicBox';

flowType.initial.type = 'uniform'; % uniform, blasius or file
flowType.initial.U0 = 0;

flowType.initial.addNoise = 1e-2;

e0 = 1/((flowParameters.gamma^2-flowParameters.gamma)*flowParameters.Ma^2);

flowType.disturb{1}.x = [-inf inf];
flowType.disturb{1}.y = [-inf inf];
flowType.disturb{1}.z = [-inf inf];
flowType.disturb{1}.var = 'UVWE';
flowType.disturb{1}.type = 'turbulenceGenerator';
flowType.disturb{1}.extraNodes = [0 0 0 0 0 0];
flowType.disturb{1}.par = [epsilon, e0];
flowType.disturb{1}.active = true;
flowType.disturb{1}.fitPoints = false;
flowType.disturb{1}.forcing = true;

%% Mesh parameters
mesh.x.n = nx;
mesh.y.n = ny;
mesh.z.n = nz;

% Available mesh types and their specific parameters are:
%   uniform
%
%   attractors
%       attractorPoints
%       attractorStrength
%       attractorSize
%
%   file
%       file


mesh.x.type = 'uniform';
mesh.x.matchFixed = 2;
mesh.x.periodic = true;
mesh.x.fixPeriodicDomainSize = true;
mesh.x.extraRefinement = 0;

mesh.y.type = 'uniform';
mesh.y.matchFixed = 2;
mesh.y.periodic = true;
mesh.y.fixPeriodicDomainSize = true;
mesh.y.extraRefinement = 0;

mesh.z.type = 'uniform';
mesh.z.matchFixed = 2;
mesh.z.periodic = true;
mesh.z.fixPeriodicDomainSize = true;
mesh.z.extraRefinement = 0;

mesh.x.buffer.i.n = 0;
mesh.x.buffer.f.n = 0;

mesh.y.buffer.i.n = 0;
mesh.y.buffer.f.n = 0;

mesh.z.buffer.i.n = 0;
mesh.z.buffer.f.n = 0;

%% Tracked points
% These are points that will have their variables written to the log file
% The mesh will also be fitted to them
%mesh.trackedPoints = [2.226810857889 0 0];

%% Time control
% If time.control = dt, qtimes and tmax are in number of iterations and the step size is fixed to dt
% If time.control = cfl, qtimes and tmax are in non-dimensional time and dt defines the maximum step size
time.control = 'cfl';

time.dt = 1;
time.maxCFL = 1.3;

time.qtimes = 1;
time.tmax = 100;

logAll = 1;

mesh.trackedPoints = [0,0,0];

mesh.trackedNorm = true;

%% Numerical methods
numMethods.spatialDerivs = 'SL6'; % SL4 or EX2
numMethods.spatialDerivsBuffer = 'EX4';
%numMethods.metricMethod = ''; % Set the method for cimputing the metrics (optional, default = 'SL4', leave blank for same as spatialDerivs)
numMethods.timeStepping = 'RK4'; % RK4 or Euler
numMethods.neumannOrder = 6;
numMethods.neumann2Order = 2;
numMethods.spatialFilterStrength = 0.499; % -0.5 < alpha < 0.5
numMethods.spatialFilterTime = 0; % Characteristic time of the spatial filter (optional, default = 0)
numMethods.filterDirections = [1 1 1];
numMethods.filterBorders = 'reducedOrder';

numMethods.SFD.type = 0; % 0 = off, 1 = whole domain, 2 = buffer zone only;
numMethods.SFD.X = 0.05;
numMethods.SFD.Delta = 10;

numMethods.SFD.applyY = false;

%numMethods.SFD.extraRegion{1}.location = [(x1+x2)/2 0 0];
%numMethods.SFD.extraRegion{1}.size = [L D inf];
%numMethods.SFD.extraRegion{1}.X = 0.1;
