function flow = generateInitialFlow(mesh,flowParameters,initialFlow,walls,flowName)
% This function defines each type of initial flow. New types may be added here if needed

nx = mesh.nx;
ny = mesh.ny;
nz = mesh.nz;

X = mesh.X;
Y = mesh.Y;
Z = mesh.Z;

gamma = flowParameters.gamma;
Ma = flowParameters.Ma;
Re = flowParameters.Re;

E0 = 1/((gamma^2-gamma)*Ma^2);

switch initialFlow.type
    case 'uniform'
        
        U = ones(nx,ny,nz);
        V = zeros(nx,ny,nz);
        W = zeros(nx,ny,nz);
        R = ones(nx,ny,nz);
        E = ones(nx,ny,nz) * E0;
        
        [~,y0ind] = min(abs(Y));
        
        if nargin == 4 || contains(flowName, 'boundaryLayer');
            U(:,1:y0ind,:) = 0;
        end

        if isfield(initialFlow,'U0')
            U = initialFlow.U0*U;
        end
        
    case 'poiseuille'
        
        U = zeros(nx,ny,nz);
        V = zeros(nx,ny,nz);
        W = zeros(nx,ny,nz);
        R = ones(nx,ny,nz);
        E = ones(nx,ny,nz) * E0;

        eta = (Y-Y(1))/(Y(end)-Y(1));

        u0 = flowParameters.U0;
        u1 = flowParameters.lowerWallVelocity;
        u2 = flowParameters.upperWallVelocity;

        U = U + (-6*u0 +3*u1 +3*u2)*eta.^2 + (6*u0 -4*u1 -2*u2)*eta + u1;

    case 'blasius'

        ybl = 0:0.0001:10;
        ubl = blasius(ybl);
        thetabl = 0.664155332943009;

        e0 = 1/((gamma^2-gamma)*Ma^2);

        nx = length(X);
        ny = length(Y);
        nz = length(Z);

        R = ones(nx,ny,nz);
        U = ones(nx,ny);
        V = zeros(nx,ny,nz);
        W = zeros(nx,ny,nz);
        E = ones(nx,ny,nz) * e0;

        % If blasiusFit defined, the blasius boundary layer will fit
        % to the wall.
        % The value of blasiusFit defines the maximum slope of Y0 for
        % the blasius profile

        Y0 = zeros(1,nx);
        if isfield(initialFlow,'blasiusFit')
            walls2D = any(walls,3);
            for i = 1:nx
                Y0(i) = Y(find(walls2D(i,:)==0,1,'first'));
            end
            for i = 2:nx
                dx = X(i)-X(i-1);
                if Y0(i-1)-Y0(i) > dx*initialFlow.blasiusFit
                    Y0(i) = Y0(i-1) - dx*initialFlow.blasiusFit;
                end
            end
            for i = nx-1:-1:1
                dx = X(i+1)-X(i);
                if Y0(i+1)-Y0(i) > dx*initialFlow.blasiusFit
                    Y0(i) = Y0(i+1) - dx*initialFlow.blasiusFit;
                end
            end
        end

        for i = 1:nx
            xi = X(i);
            if xi > 0
                theta = 0.664 * sqrt(xi/Re);
                U(i,:) = interp1(ybl*theta/thetabl+Y0(i),ubl,Y(1:end),'spline')';
                U(i,Y<Y0(i)) = 0;
            end
        end
        
        U(U>1) = 1;
        U(isnan(U)) = 1;

        U = repmat(U,[1 1 nz]);
        
    case 'compressibleBL_isothermal'
        compressibleBL_flow = calcCompressibleBL(flowParameters,false,mesh);

        U = compressibleBL_flow.U;
        V = compressibleBL_flow.V;
        W = compressibleBL_flow.W;
        E = compressibleBL_flow.E;
        R = compressibleBL_flow.R;

    case 'compressibleBL_adiabatic'
        compressibleBL_flow = calcCompressibleBL(flowParameters,true,mesh);

        U = compressibleBL_flow.U;
        V = compressibleBL_flow.V;
        W = compressibleBL_flow.W;
        E = compressibleBL_flow.E;
        R = compressibleBL_flow.R;

    case 'file'
		
		% If a folder is given, find the last flow file inside it
		if strcmp(initialFlow.flowFile(end),'/')
			nStep = checkPreviousRun(initialFlow.flowFile(1:end-1));
            if ~isempty(nStep)
                initialFlow.flowFile = sprintf('%sflow_%.10d.mat',initialFlow.flowFile,nStep);
            else
                initialFlow.flowFile = sprintf('%sbaseflow.mat',initialFlow.flowFile);
            end
			clear nStep
		end
		
		flowFile = load(initialFlow.flowFile);
		
        % If a mesh file is given, interpolate from it
        if isfield(initialFlow,'meshFile')

			% If a folder is given, use the mesh file inside it
			if strcmp(initialFlow.meshFile(end),'/')
				initialFlow.meshFile = [initialFlow.meshFile 'mesh.mat'];
			end

            meshFile = load(initialFlow.meshFile);
            Xfile = meshFile.X;
            Yfile = meshFile.Y;
            Zfile = meshFile.Z;
            
            Ufile = flowFile.U;
            Vfile = flowFile.V;
            Wfile = flowFile.W;
            Rfile = flowFile.R;
            Efile = flowFile.E;
            
            % Change NaNs to default values
            Ufile(isnan(Ufile)) = 0;
            Vfile(isnan(Vfile)) = 0;
            Wfile(isnan(Wfile)) = 0;
            Rfile(isnan(Rfile)) = 1;
            Efile(isnan(Efile)) = E0;
            
            % Create a temporary mesh and crop it to avoid extrapolations
            [Xmesh,Ymesh,Zmesh] = ndgrid(X,Y,Z);
            Xmesh(Xmesh < Xfile(1)) = Xfile(1);
            Xmesh(Xmesh > Xfile(end)) = Xfile(end);
            Ymesh(Ymesh < Yfile(1)) = Yfile(1);
            Ymesh(Ymesh > Yfile(end)) = Yfile(end);
            Zmesh(Zmesh < Zfile(1)) = Zfile(1);
            Zmesh(Zmesh > Zfile(end)) = Zfile(end);
            
            % Interpolate
            if nz == 1 || length(Zfile) == 1 % In a 2D case
                U = interpn(Xfile,Yfile,Ufile,Xmesh,Ymesh,'linear');
                V = interpn(Xfile,Yfile,Vfile,Xmesh,Ymesh,'linear');
                W = interpn(Xfile,Yfile,Wfile,Xmesh,Ymesh,'linear');
                R = interpn(Xfile,Yfile,Rfile,Xmesh,Ymesh,'linear');
                E = interpn(Xfile,Yfile,Efile,Xmesh,Ymesh,'linear');
            else % In a 3D case
                U = interpn(Xfile,Yfile,Zfile,Ufile,Xmesh,Ymesh,Zmesh,'linear');
                V = interpn(Xfile,Yfile,Zfile,Vfile,Xmesh,Ymesh,Zmesh,'linear');
                W = interpn(Xfile,Yfile,Zfile,Wfile,Xmesh,Ymesh,Zmesh,'linear');
                R = interpn(Xfile,Yfile,Zfile,Rfile,Xmesh,Ymesh,Zmesh,'linear');
                E = interpn(Xfile,Yfile,Zfile,Efile,Xmesh,Ymesh,Zmesh,'linear');
            end
			
			% Change Mach number if needed
			if isfield(initialFlow,'changeMach') && initialFlow.changeMach
				if meshFile.flowParameters.Ma ~= Ma
					fprintf('Initial flow file has a Mach number of %f which will be changed to %f to match the current simulation\n',meshFile.flowParameters.Ma,Ma);
					E = E * meshFile.flowParameters.Ma^2/Ma^2;
				end
			end

        else
            % If no mesh file is given, just load the values and check the size
            U = flowFile.U;
            if nx ~= size(U,1) || ny ~= size(U,2) || (nz ~= size(U,3) && size(U,3) ~= 1)
                error('Mesh size in initial flow file is not consistent with current mesh (%dx%dx%d -> %dx%dx%d)\nConsider changing the parameters file or providing a mesh file',size(U,1),size(U,2),size(U,3),nx,ny,nz)
            end
            V = flowFile.V;
            W = flowFile.W;
            R = flowFile.R;
            E = flowFile.E;
            
            % If the mesh provided was 2D and the case is 3D, replicate it
            if nz == 1 && size(U,3)
                U = repmat(U,[1 1 nz]);
                V = repmat(V,[1 1 nz]);
                W = repmat(W,[1 1 nz]);
                R = repmat(R,[1 1 nz]);
                E = repmat(E,[1 1 nz]);
            end
            
        end
        
end

if isfield(initialFlow,'addNoise')

	if ~isfield(initialFlow,'noiseType') || strcmp(initialFlow.noiseType,'rand')
		noiseU = initialFlow.addNoise*randn([nx ny nz]);
		noiseV = initialFlow.addNoise*randn([nx ny nz]);
		noiseW = initialFlow.addNoise*randn([nx ny nz]);
		noiseR = initialFlow.addNoise*randn([nx ny nz]);
		noiseE = initialFlow.addNoise*randn([nx ny nz]);
	elseif strcmp(initialFlow.noiseType,'uniform')
		noiseU = initialFlow.addNoise*ones([nx ny nz]);
		noiseV = initialFlow.addNoise*ones([nx ny nz]);
		noiseW = initialFlow.addNoise*ones([nx ny nz]);
		noiseR = initialFlow.addNoise*ones([nx ny nz]);
		noiseE = initialFlow.addNoise*ones([nx ny nz]);
	end
	
	if isfield(initialFlow,'noiseCenter')
		if ~isfield(initialFlow,'noiseSigma')
			initialFlow.noiseSigma = [(X(end)-X(1)).^2/10 (Y(end)-Y(1)).^2/10 (Z(end)-Z(1)).^2/10];
		end
		x0 = initialFlow.noiseCenter(1);
		y0 = initialFlow.noiseCenter(2);
		z0 = initialFlow.noiseCenter(3);
		
		sigmaX = initialFlow.noiseSigma(1);
		sigmaY = initialFlow.noiseSigma(2);
		sigmaZ = initialFlow.noiseSigma(3);
    else
        x0 = 0;
        y0 = 0;
        z0 = 0;
        sigmaX = inf;
        sigmaY = inf;
        sigmaZ = inf;
    end
    
    if nz == 1
        sigmaZ = inf;
    end
    radius = bsxfun(@plus,(X'-x0).^2/sigmaX,(Y-y0).^2/sigmaY);
    radius = bsxfun(@plus,radius,(permute(Z,[1 3 2])-z0).^2/sigmaZ);
    noiseGaussian = exp(-radius);
    
	U = U + noiseU.*noiseGaussian;
	V = V + noiseV.*noiseGaussian;
	if nz > 1
		W = W + noiseW.*noiseGaussian;
	end
	R = R + noiseR.*noiseGaussian;
	E = E + noiseE.*noiseGaussian;
end

flow.U = U;
flow.V = V;
flow.W = W;
flow.R = R;
flow.E = E;

end

function [u] = blasius(y)

f0 = [0 0 0.33204312];

[~,u] = ode45(@blasiusEq,y,f0);

u = u(:,2);

end

function d = blasiusEq(~,f)
    d(3) = -0.5*f(1)*f(3);
    d(2) = f(3);
    d(1) = f(2);
    d = d';
end
