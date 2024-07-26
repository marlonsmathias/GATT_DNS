Re = 2000;
Ma = 0.5;
L = 4;
H = 1;
dpdx = 16/Re;

%% Case name
caseName = ['Poiseuille-Re' num2str(Re) '-Ma' strrep(num2str(Ma),'.','') '-L' num2str(L) '-H' num2str(H)];

%% Domain decomposition
p_row = 10;
p_col = 1;

%% Flow parameters
flowParameters.Re = Re;
flowParameters.Ma = Ma;
flowParameters.Pr = 0.71;
flowParameters.gamma = 1.4;
flowParameters.T0 = 300;

%% Domain parameters
domain.xi = 0;
domain.xf = L;
domain.yi = -H/2;
domain.yf = H/2;
domain.zi = 0;
domain.zf = 1;

%% Flow type
flowType.name = 'poiseuilleFlow';
flowParameters.U0 = 1; % Initial mean velocity
flowParameters.lowerWallVelocity = 0;
flowParameters.upperWallVelocity = 0;

flowType.initial.type = 'uniform'; % uniform, blasius or file
flowType.initial.addNoise = 1e-4;

flowType.disturb{1}.x = [-inf inf];
flowType.disturb{1}.y = [-inf inf];
flowType.disturb{1}.z = [-inf inf];
flowType.disturb{1}.var = 'RU';
flowType.disturb{1}.type = 'pressureGradient';
flowType.disturb{1}.extraNodes = [0 0 0 0 0 0];
flowType.disturb{1}.par = [dpdx];
flowType.disturb{1}.active = true;
flowType.disturb{1}.fitPoints = false;
flowType.disturb{1}.forcing = true;

%% Mesh parameters
mesh.x.n = 100;
mesh.y.n = 200;
mesh.z.n = 1;

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
mesh.x.fixPeriodicDomainSize = false;
mesh.x.extraRefinement = 0;

mesh.y.type = 'tanh';
mesh.y.local = 'b';
mesh.y.par = 2;
mesh.y.matchFixed = 2;
mesh.y.periodic = false;
mesh.y.fixPeriodicDomainSize = false;
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

%% Time control
% If time.control = dt, qtimes and tmax are in number of iterations and the step size is fixed to dt
% If time.control = cfl, qtimes and tmax are in non-dimensional time and dt defines the maximum step size
time.control = 'cfl';

time.dt = 1;
time.maxCFL = 1.3;

time.qtimes = 1;
time.tmax = 1000;

logAll = 25;

mesh.trackedPoints = [0,0,0];

mesh.trackedNorm = true;

%% Numerical methods
numMethods.spatialDerivs = 'SL6'; % SL4 or EX2
numMethods.spatialDerivsBuffer = 'EX4';
%numMethods.metricMethod = ''; % Set the method for cimputing the metrics (optional, default = 'SL4', leave blank for same as spatialDerivs)
numMethods.timeStepping = 'RK4'; % RK4 or Euler
numMethods.neumannOrder = 6;
numMethods.neumann2Order = 2;
numMethods.spatialFilterStrength = 0.49; % -0.5 < alpha < 0.5
numMethods.spatialFilterTime = 0; % Characteristic time of the spatial filter (optional, default = 0)
numMethods.filterDirections = [1 1 1];
numMethods.filterBorders = 'reducedOrder';
numMethods.filterBordersStartX = false;
numMethods.filterBordersEndX = false;
numMethods.filterBordersEndY = false;

numMethods.SFD.type = 0; % 0 = off, 1 = whole domain, 2 = buffer zone only;
numMethods.SFD.X = 0.05;
numMethods.SFD.Delta = 10;

numMethods.SFD.applyY = false;

%numMethods.SFD.extraRegion{1}.location = [(x1+x2)/2 0 0];
%numMethods.SFD.extraRegion{1}.size = [L D inf];
%numMethods.SFD.extraRegion{1}.X = 0.1;
