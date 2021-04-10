# GATT_DNS

This is a Direct Numerical Simulation code for the compressible Navier-Stokes equations. The linear stability analysis (LST) is done by a Jacobian-free time-stepping algorithm.

The main features are listed below:

* Structured mesh (2D or 3D)
* Mesh stretching to concentrate nodes in the regions of interest
* 4th order Runge-Kutta for the time-stepping
* 4th order spectral-like spatial differentiation
* 10th order spatial anti-aliasing filter
* Temporal low-pass filter (SFD) to allow reaching a base-flow even ate unstable conditions
* Buffer-zones at the open boundaries
* Domain decomposition for parallel execution

The code was developed by the Group of Aeroacoustics, Transition and Turbulence of the São Carlos School of Engineering - University of São Paulo (EESC-USP)

To cite this work, please use the following paper, which brings more details on the implementation and some results obtained by the code:  
Mathias MS, Medeiros M. Direct Numerical Simulation of a Compressible Flow and Matrix-Free Analysis of its Instabilities over an Open Cavity. J Aerosp Technol Manag [Internet]. 2018 Jul 20;10:1–13. Available from: http://www.jatm.com.br/ojs/index.php/jatm/article/view/949

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# System requirements:

Matlab (tested on 2018a)

gfortran (tested on 9.3.0)

openmpi (tested on 4.0.3, libopenmpi-dev package)

2decomp&fft (tested on 1.5, compiled with the same compiler as the code, http://www.2decomp.org/)


# Main files:

runDNS.m

This file is the one that will call all the others. The parameters file name should be defined in it.

The first output contains a handle to each of the saved flows.

The second output contains all internal variables created by Matlab at the runtime.

The parameters to be chosen are:

caseFile -> Name of the parameters file to be called. Do not include the .m extension.

runSimulation -> Whether or the the DNS will be run

plotDNSDomain -> Whether or not the domain will be plotted before runtime

forceRecompile -> Force the recompilation to be triggered, even if nothing has changed from previous runs.

forceRecompileAll -> Force all files to be recompiled instead of just the ones that changed.

displayCompiling -> Whether or not the compiler output is shown in the screen.

optimizeCode -> Whether or not optimization flags will be passed to the compiler. If true, compilation will take longer but the code will run much faster.

debugger -> Whether or not the gdb debugger is called at runtime. This overrides the optimization flag. The X server must be on, as one xterm window will open for each process.

profiler -> Create a profile.txt file with code performance data. gprof is used. This increases the total runtime.

matlabDir -> Root directory of the Matlab installation. Leave empty for auto.

decompDir -> Directory where the 2decomp library is found.


## Parameters file

parameters.m (The actual name is chosen in runDNS.m or passed as an argument)

This file contains all the parameters for the simulation, both for the flow and for the numerical methods.

A copy of this file is made at caseName/Fortran at runtime.

## Example parameters file:

	%% Case name
	caseName = 'testRun'; %This is the name of the folder the results are saved to.

	%% Domain decomposition
	p_row = 4;
	p_col = 1;

	%% Flow parameters
	flowParameters.Re = 1000;
	flowParameters.Ma = 0.5;
	flowParameters.Pr = 0.71;
	flowParameters.gamma = 1.4;
	flowParameters.T0 = 300;

	%% Domain parameters
	domain.xi = -1;
	domain.xf = 5;
	domain.yi = -1;
	domain.yf = 4;
	domain.zi = 0;
	domain.zf = 1;

	%% Flow type
	flowType.name = 'boundaryLayerIsothermal'; % Check the source/boundaries folder for available flow types.

	flowType.initial.type = 'uniform'; % Either uniform, blasius or file. This is the initial flow that is created if no previous run is found.

	%flowType.initial.blasiusFit = 1; % For initial type of blasius. Fits the bottom of the profile to the geometry, with a maximum slope of blasiusFit (optional)

	%flowType.initial.flowFile = 'baseflows/TSWaves/flow.mat';
	%flowType.initial.meshFile = 'baseflows/TSWaves/mesh.mat'; % Providing a mesh is optional. If no mesh is given the code assumes it is the same as the current mesh. If the given mesh is different, the flow will be interpolated to the new mesh.
	%flowType.initial.addNoise = 1e-5; % Adds white noise of this magnitude to the initial flow (optional)
	%flowType.initial.changeMach = true; % Changes the Mach number in the initial flow file if it does not match the current Mach number (optional)


	% Positions of the cavities
	flowType.cav{1}.x = [1 2];
	flowType.cav{1}.y = [-1 0];
	flowType.cav{1}.z = [-inf inf];

	% Positions of the roughnesses
	%flowType.rug{1}.x = [1 1.2];
	%flowType.rug{1}.y = [-0.1 0];
	%flowType.rug{1}.z = [0.3 0.8];
	 
	%flowType.rug{2}.x = [2.5 4];
	%flowType.rug{2}.y = [0.5 1];
	%flowType.rug{2}.z = [0.3 0.8];

	% Positions of the disturbances
	%flowType.disturb{1}.x = [0.1 0.2];
	%flowType.disturb{1}.y = [0 0];
	%flowType.disturb{1}.z = [-inf inf];
	%flowType.disturb{1}.var = 'V'; % Which variables it applies to, among UVWRE.
	%flowType.disturb{1}.type = 'periodicSine'; % Name of the disturbances file. Check the source/disturbances folder.
	%flowType.disturb{1}.extraNodes = [0 0 0 0 0 0]; % Amount of extra nodes that will be passed to the subroutine, in the following order: xi xf yi yf zi zf (optional)
	%flowType.disturb{1}.par = [1, 1e-4]; % Extra parameters to the passed to the disturbance routine. (optional)
	%flowType.disturb{1}.active = true;
	%flowType.disturb{1}.fitPoints = true; % Whether or not the mesh will be fitted to this disturbance. (optional)

	%% Mesh parameters
	mesh.x.n = 300;
	mesh.y.n = 200;
	mesh.z.n = 1;

	%mesh.x.d0 = 0.01; % If mesh.n is not provided, mesh.d0 will be used to compute n so that the base spacing of the mesh is d0

	% Available mesh types and their specific parameters are:
	%   uniform
	%
	%   power
	%       power
	%
	%   attractors
	%       attractorPoints
	%       attractorStrength
	%       attractorSize
			attractorRegions (optional)
	%
	%   tanh
	%       local -> i, f, b
	%       par
	%
	%   file
	%       file
	%       fileCalcBuffer (optional)


	%mesh.x.type = 'uniform';
	mesh.x.type = 'attractors';
	mesh.x.attractorPoints = [0 1 2];
	mesh.x.attractorStrength = [1 1 2];
	mesh.x.attractorSize = [0.3 1 1];
	%mesh.x.attractorRegions = [3 4 5 6 1]; % Mesh refines from x1 to x2, then is constant until x3 and coarses back until x4. x5 is the strength. Each line describes one region.
	%mesh.x.type = 'file';
	%mesh.x.file = 'meshes/x.dat';
	%mesh.x.fileCalcBuffer = true; % (optional, default = false) Set to true if the mesh file does not contain the buffer zone
	mesh.x.matchFixed = true; % If true, the mesh may be slightly transformed so that nodes are present exactly at geometry corners and flow disturbances. If matchFixed=2, a pchip interpolation will be used instead of spline.
	mesh.x.periodic = false;
	mesh.x.fixPeriodicDomainSize = false; % If true, the periodic domain length will be slight reduced so that the space between nodes n and 1 is accounted for.
	mesh.x.extraRefinement = 0; % This will add n nodes between each pair of nodes. Useful for mesh refinement tests.

	%mesh.y.type = 'power';
	%mesh.y.power = 2;
	mesh.y.type = 'attractors';
	mesh.y.attractorPoints = 0;
	mesh.y.attractorStrength = 10;
	mesh.y.attractorSize = 0.5;
	%mesh.y.type = 'file';
	%mesh.y.file = 'meshes/y.dat';
	mesh.y.matchFixed = true;
	mesh.y.periodic = false;
	mesh.y.fixPeriodicDomainSize = false;
	mesh.y.extraRefinement = 0;

	mesh.z.type = 'uniform';
	mesh.z.matchFixed = true;
	mesh.z.periodic = true;
	mesh.z.fixPeriodicDomainSize = true;
	mesh.z.extraRefinement = 0;

	mesh.x.buffer.i.n = 20;
	mesh.x.buffer.i.type = 'sigmoid'; % Sigmoid or exponential
	mesh.x.buffer.i.stretching = 20;

	mesh.x.buffer.f.n = 20;
	%mesh.x.buffer.f.l = 300; % Buffer zone length can be specifiend instead of number of nodes
	mesh.x.buffer.f.type = 'sigmoid';
	mesh.x.buffer.f.stretching = 20;
	mesh.x.buffer.f.transition = 0.2;

	mesh.y.buffer.i.n = 0;

	mesh.y.buffer.f.n = 20;
	mesh.y.buffer.f.type = 'exponential';
	mesh.y.buffer.f.stretching = 0.1;
	mesh.y.buffer.f.ramp = 20; % Defines number of nodes used to ramp up the spacing change in an exponential buffer zone (optional)

	mesh.z.buffer.i.n = 0;
	mesh.z.buffer.f.n = 0;

	%% Tracked points
	% These are points that will have their variables written to the log file
	% The mesh will also be fitted to them
	mesh.trackedPoints = [1 1 0; 2 1.5 0];
	mesh.fitTrackedPoints = false; % If true, the mesh will be fitted to these points. (optional)
	mesh.trackedNorm = true; % Normalize the values of probe points in the log file (optional, default = false).

	logAll = false; % If true, all step will be saved to the log file. (optional)

	%% Time control
	% If time.control = dt, qtimes and tmax are in number of iterations and the step size is fixed to dt
	% If time.control = cfl, qtimes and tmax are in non-dimensional time and dt defines the maximum step size
	time.control = 'dt';

	time.dt = 1.5e-3;
	time.maxCFL = 1;

	time.qtimes = 100;
	time.tmax = 1000;

	%time.CFLignoreZ = true % Ignores the Z direction for CFL (optional, default = false)

	%% Numerical methods
	numMethods.spatialDerivs = 'SL4'; % SL4, EX2, EX4 or EXn (with n = order)
	numMethods.spatialDerivsBuffer = 'EX2';
	numMethods.timeStepping = 'RK4'; % RK4, Euler or SSPRK3
	numMethods.neumannOrder = 4;
	numMethods.neumann2Order = 2;
	numMethods.spatialFilterStrength = 0.49; % -0.5 < alpha < 0.5
	numMethods.spatialFilterTime = 0.01; % Characteristic time of the spatial filter (optional, default = 0)
	numMethods.filterDirections = [1 1 1];
	numMethods.filterBorders = true; ('decentered' (default) or 'reducedOrder' for the strategy near the wall)
	numMethods.filterBordersStartX = false; % Whether or not to use the spatial filter at the start of the X domain (optional, default = false)
	numMethods.filterBordersEndX = false;
	numMethods.filterBordersEndY = false;

	%numMethods.changeOrderZ = false; % Keep the spatial derivative order in the buffer zone for each direction (optional, default = true)

	numMethods.SFD.type = 0; % 0 = off, 1 = whole domain, 2 = buffer zone only;
	numMethods.SFD.X = 0.5;
	numMethods.SFD.Delta = 20;

	%numMethods.SFD.applyZ = false; % Turn off SFD in certain buffer zones (optional, default = true);

	% This is an optional part that creates an extra region where the SFD will be active. (Only possible if SFD.type = 2)
	%numMethods.SFD.extraRegion{1}.location = [0 0 0];
	%numMethods.SFD.extraRegion{1}.size = [50 1 inf];
	%numMethods.SFD.extraRegion{1}.X = 30;

## Other Matlab files:

preprocessing.m -> Calls all other preprocessing files

compileFortran.m -> Calls the Fortran compliler with the correct flags

### Derivatives and filter matrices:

findDerivativeRegions.m -> Creates a map of the different regions for the derivatives

finiteDifferenceCoefficients.m -> Contains the coefficients for various methods

getMatrixTypeBlocks.m -> Transforms the region map into a list of blocks

makeMatrices.m -> Builds the matrices for derivatives and filters, also applies the metrics

prepareThomas.m -> Pre compute the Thomas algorithm coefficients from the matrices

spatialFilterCoefficients.m -> Computes the spatial filter coefficients

## Boundary conditions:

findWallsForBoundaries.m -> Finds the locations of the walls in the domain

getBoundaryConditions.m -> Calls the correct boundary condition routine from the source/boundaries folder

initBoundaries.m -> Separates the boundaries in blocks for the different Fortran processes.

### Initial flow:

checkPreviousRun.m -> Checks if there are any previous run files in the folder.

generateInitialFlow.m -> Contains the functions to generate the various types of initial flows.

initialFlow.m -> Gets the initial flow and saves to a file.

### Mesh:

generateMesh.m -> Contains the functions for the various types of mesh

meshAddFixedPoints.m -> Finds the fixed points where a node in required

### Write Fortran files:

writeFortranBoundaries.m

writeFortranDisturbances.m

writeFortranMatrices.m

writeFortranParameters.m

### Misc:

calcSFDregion.m -> Generates a mat file with the SFD coefficients for each node.

getDomainSlices.m -> Computes the domain decomposition in the same way 2decomp does.

plotDomain.m -> Shows the current domain in a plot window.

## Disturbance files:

Must be placed in the source/disturbances folder and be named as disturbanceName.F90. An example is provided below:

    subroutine disturbanceName(nx,ny,nz,x,y,z,t,u,par1,par2)

    implicit none
    
    integer, intent(in) :: nx, ny, nz
    real*8, dimension(nx), intent(in) :: x
    real*8, dimension(ny), intent(in) :: y
    real*8, dimension(nz), intent(in) :: z
    real*8, intent(in) :: t
    real*8, dimension(nx,ny,nz),intent(inout) :: u
    real*8, intent(in) :: par1
    real*8, intent(in) :: par2
    
    ! Modify u

    end subroutine

Multiple variables can be given as inputs. Inputs par1, par2 and etc are optional.


# Acknowledgements

The authors would like to thank the São Paulo Research Foundation (FAPESP/Brazil), for grants 2018/04584-0 and 2017/23622-8; and the National Council for Scientific and Technological Development (CNPq/Brazil) for grants 134722/2016-7 and 307956/2019-9.
