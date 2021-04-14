%% Case name
caseName = 'Red600-D5-Ma03-LD3';

%delta = 1;
%Red = 600;
%LD = 3;
%D = 5;
%x1 = delta^2*Red/(1.72^2*delta^2);
%x2 = x1 + LD*D;

%% Domain decomposition
p_row = 4;
p_col = 1;

%% Flow parameters
flowParameters.Re = 600;
flowParameters.Ma = 0.3;
flowParameters.Pr = 0.71;
flowParameters.gamma = 1.4;
flowParameters.T0 = 300;

%% Domain parameters
domain.xi = -50;
domain.xf = 500;
domain.yi = -5;
domain.yf = 50;
domain.zi = 0;
domain.zf = 1;

%% Flow type
flowType.name = 'boundaryLayerIsothermal';

flowType.initial.type = 'blasius'; % uniform, blasius or file
%flowType.initial.blasiusFit = 0.2;
%flowType.initial.flowFile = 'ReD3000-Ddelta5-Ma03-LD3-large/';
%flowType.initial.meshFile = 'ReD3000-Ddelta5-Ma03-LD3/';
flowType.initial.addNoise = 1e-3; % Adds noise of this magnitude to the initial flow (optional)
flowType.initial.noiseType = 'rand'; % Noise type, either 'rand' or 'uniform' (optional, default = 'rand')
flowType.initial.noiseCenter = [200 0 0]; % Center point of noise gaussian (optional)
flowType.initial.noiseSigma = [10 5 inf]; % Size of gaussian for noise (optional)

flowType.cav{1}.x = 5*[40.573788325646099 43.573788325646099];
flowType.cav{1}.y = 5*[-1 0];
flowType.cav{1}.z = [-inf inf];

%flowType.rug{1}.x = [1 1.2];
%flowType.rug{1}.y = [-0.1 0];
%flowType.rug{1}.z = [0.3 0.8];

%flowType.rug{2}.x = [2.5 4];
%flowType.rug{2}.y = [0.5 1];
%flowType.rug{2}.z = [0.3 0.8];

flowType.disturb{1}.x = 5*[5 10];
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
mesh.x.d0 = 2.5;
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
mesh.x.attractorRegions = [x1-L x1 x2 x2+2*L 10
                           x1-L/4 x1 x1+L/4 x1+L/2 5
						   x2-L/2 x2-L/4 x2 x2+L/2 10];
%mesh.x.type = 'file';
%mesh.x.file = 'baseflows/3Dtest/x.dat';
mesh.x.matchFixed = true;
mesh.x.periodic = false;
mesh.x.fixPeriodicDomainSize = false;
mesh.x.extraRefinement = 0;

%mesh.y.type = 'power';
%mesh.y.power = 2;
mesh.y.type = 'attractors';
mesh.y.attractorPoints = [0.5];
mesh.y.attractorStrength = [10];
mesh.y.attractorSize = [2];
mesh.y.attractorRegions = [-10 -5 1 10 10];
%mesh.y.type = 'file';
%mesh.y.file = 'baseflows/3Dtest/y.dat';
mesh.y.matchFixed = true;
mesh.y.periodic = false;
mesh.y.fixPeriodicDomainSize = false;
mesh.y.extraRefinement = 0;

mesh.z.type = 'uniform';
%mesh.z.file = 'baseflows/3Dtest/z.dat';
mesh.z.matchFixed = true;
mesh.z.periodic = true;
mesh.z.fixPeriodicDomainSize = true;
mesh.z.extraRefinement = 0;

%mesh.x.buffer.i.n = 15;
mesh.x.buffer.i.l = 300;
mesh.x.buffer.i.type = 'sigmoid';
mesh.x.buffer.i.stretching = 30;
mesh.x.buffer.i.transition = 0.2;

mesh.x.buffer.f.n = 25;
mesh.x.buffer.f.type = 'sigmoid';
mesh.x.buffer.f.stretching = 30;
mesh.x.buffer.f.transition = 0.2;

mesh.y.buffer.i.n = 0;

mesh.y.buffer.f.l = 300;
%mesh.y.buffer.f.type = 'sigmoid';
%mesh.y.buffer.f.stretching = 50;
mesh.y.buffer.f.transition = 0.2;
mesh.y.buffer.f.type = 'exponential';
mesh.y.buffer.f.stretching = 0.1;
mesh.y.buffer.f.ramp = 20;

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

time.qtimes = 50;
time.tmax = 10000;

logAll = 100;

mesh.trackedPoints = 5*[7.5 -0.2 0;
                      20 0 0;
					  25 0 0;
					  30 0 0;
					  35 0 0;
					  40 0 0;
					  45 0 0;
					  50 0 0;
					  55 0 0;
					  60 0 0;
					  65 0 0;
					  70 0 0;
					  75 0 0;
					  80 0 0;
					  85 0 0;
					  90 0 0;
					  95 0 0;
					  100 0 0] + [0 1 0];
					  
mesh.trackedPoints = [];

mesh.trackedNorm = true;


%% Numerical methods
numMethods.spatialDerivs = 'SL6'; % SL4 or EX2
numMethods.spatialDerivsBuffer = 'EX4';
numMethods.timeStepping = 'RK4'; % RK4 or Euler
numMethods.neumannOrder = 6;
numMethods.neumann2Order = 2;
numMethods.spatialFilterStrength = 0.49; % -0.5 < alpha < 0.5
numMethods.spatialFilterTime = 0; % Characteristic time of the spatial filter (optional, default = 0)
numMethods.filterDirections = [1 1 1];
numMethods.filterBorders = 'reducedOrder';
%numMethods.filterBordersStartX = false;
%numMethods.filterBordersEndX = false;
%numMethods.filterBordersEndY = false;

numMethods.SFD.type = 2; % 0 = off, 1 = whole domain, 2 = buffer zone only;
numMethods.SFD.X = 0.02;
numMethods.SFD.Delta = 50;

numMethods.SFD.applyY = false;

%numMethods.SFD.extraRegion{1}.location = [750 0 0];
%numMethods.SFD.extraRegion{1}.size = [500 5 inf];
%numMethods.SFD.extraRegion{1}.X = 0.02;
