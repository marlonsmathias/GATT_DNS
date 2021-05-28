delta = 1;
Red = 600;
Ma = 0.5;
L = 10;
D = 5;
x1 = delta^2*Red/(1.72^2*delta^2);
x2 = x1 + L;
xEnd = x2+200;
delta99end = 5*xEnd/sqrt(Red*xEnd);
yEnd = 5*delta99end;

%% Case name
caseName = ['Red' num2str(Red) '-Ma' strrep(num2str(Ma),'.','') '-L' num2str(L) '-D' num2str(D)];

%% Domain decomposition
p_row = 10;
p_col = 1;

%% Flow parameters
flowParameters.Re = Red;
flowParameters.Ma = Ma;
flowParameters.Pr = 0.71;
flowParameters.gamma = 1.4;
flowParameters.T0 = 300;

%% Domain parameters
domain.xi = 0;
domain.xf = xEnd;
domain.yi = -D;
domain.yf = yEnd;
domain.zi = 0;
domain.zf = 1;

%% Flow type
flowType.name = 'boundaryLayerIsothermal';

flowType.initial.type = 'blasius'; % uniform, blasius or file
%flowType.initial.blasiusFit = 0;
%flowType.initial.flowFile = 'Red600-Ma05-L10-D5-ppp13-refCav2q/';
%flowType.initial.meshFile = 'Red600-Ma05-L10-D5-ppp13-refCav2q/';

flowType.cav{1}.x = [x1 x2];
flowType.cav{1}.y = [-D 0];
flowType.cav{1}.z = [-inf inf];

%flowType.rug{1}.x = [1 1.2];
%flowType.rug{1}.y = [-0.1 0];
%flowType.rug{1}.z = [0.3 0.8];
 
%flowType.rug{2}.x = [2.5 4];
%flowType.rug{2}.y = [0.5 1];
%flowType.rug{2}.z = [0.3 0.8];

flowType.disturb{1}.x = [25 50];
flowType.disturb{1}.y = [0 0];
flowType.disturb{1}.z = [-inf inf];
flowType.disturb{1}.var = 'V';
flowType.disturb{1}.type = 'packet_2d';
flowType.disturb{1}.extraNodes = [0 0 0 0 0 0];
flowType.disturb{1}.par = [0.02, 50, 1e-5]; %omega,nmodes,amplitude
flowType.disturb{1}.active = false;
flowType.disturb{1}.fitPoints = false;

%% Mesh parameters
%mesh.x.n = 600;
%mesh.y.n = 120;
mesh.x.d0 = 4;
mesh.y.d0 = 1;
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


%mesh.x.type = 'uniform';
mesh.x.type = 'attractors';
x1 = flowType.cav{1}.x(1);
x2 = flowType.cav{1}.x(2);
L = x2 - x1;
mesh.x.attractorPoints = [];
mesh.x.attractorStrength = [];
mesh.x.attractorSize = [];
mesh.x.attractorRegions = [x1-100 x1-50 x2+50 x2+100 3;
                           x1-20 x1 x2+10 x2+20 16;
						   x2-5 x2-0.1 x2 x2+2 30
						   x1-2 x1 x1+0.1 x1+5 20
						   -30 -10 10 30 1];
%mesh.x.type = 'file';
%mesh.x.file = 'baseflows/3Dtest/x.dat';
mesh.x.matchFixed = 2;
mesh.x.periodic = false;
mesh.x.fixPeriodicDomainSize = false;
mesh.x.extraRefinement = 0;

%mesh.y.type = 'power';
%mesh.y.power = 2;
mesh.y.type = 'attractors';
mesh.y.attractorPoints = [];
mesh.y.attractorStrength = [];
mesh.y.attractorSize = [];
mesh.y.attractorRegions = [-2*D -D/10 D/10 2*delta99end 12];
%mesh.y.type = 'file';
%mesh.y.file = 'baseflows/3Dtest/y.dat';
mesh.y.matchFixed = 2;
mesh.y.periodic = false;
mesh.y.fixPeriodicDomainSize = false;
mesh.y.extraRefinement = 0;

mesh.z.type = 'uniform';
%mesh.z.file = 'baseflows/3Dtest/z.dat';
mesh.z.matchFixed = 2;
mesh.z.periodic = true;
mesh.z.fixPeriodicDomainSize = true;
mesh.z.extraRefinement = 0;

mesh.x.buffer.i.n = 40;
mesh.x.buffer.i.type = 'sigmoid';
mesh.x.buffer.i.stretching = 30;
mesh.x.buffer.i.transition = 0.2;
%mesh.x.buffer.i.ramp = 20;

mesh.x.buffer.f.n = 30;
mesh.x.buffer.f.type = 'sigmoid';
mesh.x.buffer.f.stretching = 15;
mesh.x.buffer.f.transition = 0.2;
%mesh.x.buffer.f.ramp = 20;

mesh.y.buffer.i.n = 0;

mesh.y.buffer.f.n = 40;
mesh.y.buffer.f.type = 'sigmoid';
mesh.y.buffer.f.stretching = 30;
mesh.y.buffer.f.transition = 0.2;
%mesh.y.buffer.f.ramp = 20;
%mesh.y.buffer.f.type = 'exponential';
%mesh.y.buffer.f.stretching = 0.1;

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

time.qtimes = 100;
time.tmax = 10000;

logAll = 25;

trackedX = 0:50:domain.xf;
nProbes = length(trackedX);
mesh.trackedPoints = [37.5 0 0; (x1+3*x2)/4 0 0; trackedX' ones(nProbes,1) zeros(nProbes,1)];

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

numMethods.SFD.type = 2; % 0 = off, 1 = whole domain, 2 = buffer zone only;
numMethods.SFD.X = 0.05;
numMethods.SFD.Delta = 10;

numMethods.SFD.applyY = false;

%numMethods.SFD.extraRegion{1}.location = [(x1+x2)/2 0 0];
%numMethods.SFD.extraRegion{1}.size = [L D inf];
%numMethods.SFD.extraRegion{1}.X = 0.1;
