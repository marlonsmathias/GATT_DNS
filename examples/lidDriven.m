Re = 1000;
Ma = 0.5;
L = 1;
D = 1;

%% Case name
caseName = ['LidDriven-Re' num2str(Re) '-Ma' strrep(num2str(Ma),'.','') '-L' num2str(L) '-D' num2str(D)];

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
domain.yi = 0;
domain.yf = D;
domain.zi = 0;
domain.zf = 1;

%% Flow type
flowType.name = 'lidDrivenFlow';

flowType.initial.type = 'uniform'; % uniform, blasius or file
flowType.initial.U0 = 0;

%% Mesh parameters
mesh.x.n = 100;
mesh.y.n = 100;
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

mesh.x.type = 'tanh';
mesh.x.local = 'b';
mesh.x.par = 2;
mesh.x.matchFixed = 2;
mesh.x.periodic = false;
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

% Add lid driven movement
flowType.disturb{1}.x = [0 L];
flowType.disturb{1}.y = [D D];
flowType.disturb{1}.z = [-inf inf];
flowType.disturb{1}.var = 'U';
flowType.disturb{1}.type = 'lidDrivenMovement';
flowType.disturb{1}.extraNodes = [0 0 0 0 0 0];
flowType.disturb{1}.par = [1]; %U0
flowType.disturb{1}.active = true;
flowType.disturb{1}.fitPoints = true;

%% Time control
% If time.control = dt, qtimes and tmax are in number of iterations and the step size is fixed to dt
% If time.control = cfl, qtimes and tmax are in non-dimensional time and dt defines the maximum step size
time.control = 'cfl';

time.dt = 1;
time.maxCFL = 1.3;

time.qtimes = 1;
time.tmax = 1000;

logAll = 25;

mesh.trackedPoints = [];

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
