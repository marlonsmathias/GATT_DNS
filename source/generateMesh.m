function [X, mesh, nx] = generateMesh(xi,xf,mesh,direction)

% This function computes the mesh for each direction and adds the buffer zones.

% If the number of nodes is undefined, generate abitrary mesh with same refining and define it
if ~isfield(mesh,'n')
    meshTemp = mesh;
	meshTemp.n = ceil((xf-xi)/mesh.d0);
    meshTemp.matchFixed = 0;
    meshTemp.extraRefinement = 0;
    meshTemp.fixPeriodicDomainSize = 0;
    meshTemp.buffer.i.n = 0;
    meshTemp.buffer.f.n = 0;
	[Xtemp] = generateMesh(xi,xf,meshTemp,direction);
    
    dBase = max(diff(Xtemp));
    
    mesh.n = ceil(meshTemp.n*dBase/mesh.d0);
end

% If the number of nodes in any buffer zone is undefined, define it based on the prescribed length
if ~isfield(mesh.buffer.i,'n') || ~isfield(mesh.buffer.f,'n')
    meshTemp = mesh;
	meshTemp.matchFixed = 0;
    meshTemp.extraRefinement = 0;
    meshTemp.fixPeriodicDomainSize = 0;
    meshTemp.buffer.i.n = 0;
    meshTemp.buffer.f.n = 0;
	
	[Xtemp] = generateMesh(xi,xf,meshTemp,direction);
	
    if ~isfield(mesh.buffer.i,'n')
        dBase = Xtemp(2) - Xtemp(1);
        target = mesh.buffer.i.l;
        mesh.buffer.i.n = floor(target/dBase);
        XB = calcBufferZone(mesh.buffer.i);
        value = XB(end)*dBase;
        
        while value > target && mesh.buffer.i.n > 2
            mesh.buffer.i.n = floor(mesh.buffer.i.n/2);
            XB = calcBufferZone(mesh.buffer.i);
            value = XB(end)*dBase;
        end
        while value < target
            mesh.buffer.i.n = mesh.buffer.i.n + 1;
            XB = calcBufferZone(mesh.buffer.i);
            value = XB(end)*dBase;
        end
    end
    if ~isfield(mesh.buffer.f,'n')
        dBase = Xtemp(end) - Xtemp(end-1);
        target = mesh.buffer.f.l;
        mesh.buffer.f.n = floor(target/dBase);
        XB = calcBufferZone(mesh.buffer.f);
        value = XB(end)*dBase;
        
        while value > target && mesh.buffer.f.n > 2
            mesh.buffer.f.n = floor(mesh.buffer.f.n/2);
            XB = calcBufferZone(mesh.buffer.f);
            value = XB(end)*dBase;
        end
        while value < target
            mesh.buffer.f.n = mesh.buffer.f.n + 1;
            XB = calcBufferZone(mesh.buffer.f);
            value = XB(end)*dBase;
        end
    end
    
end

% If mesh is a single point, run as uniform
if mesh.n == 1
	mesh.type = 'uniform';
end

% Remove buffer zone for periodic cases
if mesh.periodic && (mesh.buffer.i.n > 0 || mesh.buffer.f.n > 0)
    warning(['Buffer zone was removed for periodic dimension ' direction]);
    mesh.buffer.i.n = 0;
    mesh.buffer.f.n = 0;
end

% Add temporary node at the end if needed
if mesh.periodic && mesh.fixPeriodicDomainSize && mesh.n > 1 && ~strcmp(mesh.type,'file')
    mesh.n = mesh.n+1;
    addedTempNode = true;
else
    addedTempNode = false;
end

% Count nodes
nx = mesh.n + mesh.buffer.i.n + mesh.buffer.f.n;
physicalStart = mesh.buffer.i.n + 1;
physicalEnd = mesh.buffer.i.n + mesh.n;

%% Compute physical domain
switch mesh.type
    case 'uniform'
        XPhysical = linspace(xi,xf,mesh.n);

    case 'power'
        XPhysical = xi + (xf-xi)*linspace(0,1,mesh.n).^mesh.power;

    case 'tanh'
        switch mesh.local
            case 'f' % Refine at the end
                eta = tanh(linspace(0,mesh.par,mesh.n));
                eta = eta/eta(end);
                
            case 'i' % Refine at the start
                eta = tanh(linspace(0,mesh.par,mesh.n));
                eta = eta/eta(end);
                eta = 1-eta(end:-1:1);
            
            case 'b' % Refine at both sides
                eta = tanh(linspace(-mesh.par,mesh.par,mesh.n));
                eta = eta/eta(end);
                eta = (eta+1)/2;
        end
        
        XPhysical = (xf-xi)*eta+xi;
        
    case 'attractors'
        xBase = linspace(xi,xf,100*mesh.n);
        eta = ones(1,100*mesh.n);
        for i = 1:length(mesh.attractorPoints)
            eta = eta + mesh.attractorStrength(i)*exp(-((xBase-mesh.attractorPoints(i))/mesh.attractorSize(i)).^2);
        end
        if isfield(mesh,'attractorRegions')
            for i = 1:size(mesh.attractorRegions,1)
                nodePositions = mesh.attractorRegions(i,1:4);
                if isinf(nodePositions(2))
                    nodePositions(1:2) = [xi-2 xi-1];
                end
                if isinf(nodePositions(3))
                    nodePositions(3:4) = [xf+1 xf+2];
                end
                nodePositions = [nodePositions(1)-1 nodePositions nodePositions(4)+1];
                eta = eta + mesh.attractorRegions(i,5)*interp1(nodePositions,[0 0 1 1 0 0],xBase,'pchip').^2;
            end
        end
        eta = cumsum(eta);
        eta = eta - eta(1);
        eta = (mesh.n-1)*eta/eta(end) + 1;
        eta(end) = mesh.n;

        XPhysical = interp1(eta,xBase,1:mesh.n,'spline');
		
    case 'attractors_old'
        xBase = linspace(xi,xf,mesh.n);
        eta = ones(1,mesh.n);
        for i = 1:length(mesh.attractorPoints)
            eta = eta + mesh.attractorStrength(i)*exp(-((xBase-mesh.attractorPoints(i))/mesh.attractorSize(i)).^2);
        end
        eta = cumsum(eta);
        eta = eta - eta(1);
        eta = (mesh.n-1)*eta/eta(end) + 1;
        eta(end) = mesh.n;

        XPhysical = interp1(eta,xBase,1:mesh.n,'spline');
        
    case 'file'
        
        X = load(mesh.file);
        if strcmp(mesh.file(end-3:end),'.mat')
            eval(['X = X.' direction ';'])
        end
        if size(X,2) == 1
            X = X';
        end
        
        if ~isfield(mesh,'fileCalcBuffer') || ~mesh.fileCalcBuffer
            if length(X) ~= nx
                error('%s mesh from file %s contains %d nodes instead of %d as specified in parameters',direction,mesh.file,length(X),nx)
            end
            return
        end
        
        if length(X) ~= mesh.n
            error('%s mesh from file %s contains %d nodes instead of %d as specified in parameters',direction,mesh.file,length(X),mesh.n)
        end
        XPhysical = X;
end
    
%% Match fixed points such as cavity edges
if mesh.matchFixed && mesh.n > 1
    fixPoints = unique([XPhysical([1 2]) mesh.fixPoints XPhysical([end-1 end])]);
    [~,closestNodes] = min(abs(bsxfun(@minus,XPhysical',fixPoints)),[],1);

    ok = ~any(diff(closestNodes) == 0); % Check for two fixed points assigned to the same node
    iter = 0;
    iterMax = mesh.n;
    while ~ok && iter <= iterMax
        duplicates = find(diff(closestNodes) == 0);
        closestNodes(duplicates) = closestNodes(duplicates) - 1;
        ok = ~any(diff(closestNodes) == 0);
        iter = iter + 1;
    end

	if mesh.matchFixed == 2
		XPhysical = XPhysical - interp1(closestNodes,XPhysical(closestNodes)-fixPoints,1:mesh.n,'pchip');
	else
		XPhysical = XPhysical - interp1(closestNodes,XPhysical(closestNodes)-fixPoints,1:mesh.n,'spline');
	end
    
    XPhysical(closestNodes) = fixPoints; % Done to fix any possible rounding errors
end

X = zeros(1,nx);
X(physicalStart:physicalEnd) = XPhysical;

%% Add buffer zones

if mesh.buffer.i.n > 0
    baseDist = X(physicalStart+1) - X(physicalStart);
    XB = calcBufferZone(mesh.buffer.i) * baseDist;
    
    X(1:physicalStart-1) = X(physicalStart) - XB(end:-1:1);
end

if mesh.buffer.f.n > 0
    baseDist = X(physicalEnd) - X(physicalEnd-1);
    XB = calcBufferZone(mesh.buffer.f) * baseDist;
    
    X(physicalEnd+1:end) = X(physicalEnd) + XB;
end

%% Add extra refinement nodes
if mesh.extraRefinement > 0
    
    % Recalculate number of nodes
    er = mesh.extraRefinement;
    mesh.n = (1+er)*mesh.n-er;
    mesh.buffer.i.n = (1+er)*mesh.buffer.i.n;
    mesh.buffer.f.n = (1+er)*mesh.buffer.f.n;
    nx = mesh.n + mesh.buffer.i.n + mesh.buffer.f.n;
    
    % Add new nodes
    Xnew = interp1(linspace(0,1,length(X)),X,linspace(0,1,nx),'pchip');
    Xnew(1:(1+er):nx) = X;
    X = Xnew;
end

%% Remove temporary node
if addedTempNode
    X(end) = [];
    mesh.n = mesh.n-1;
    nx = nx-1;
end

end

function XB = calcBufferZone(par)

    switch par.type
        case 'exponential'
        
			if ~isfield(par,'ramp')
				XB(1) = 1+par.stretching;
				for i = 2:par.n
					XB(i) = XB(i-1) + (1+par.stretching)^i;
				end
			
			else
				delta = ones(1,par.n);
				for i = 2:par.n
					stretching = 1 + min(1,i/par.ramp) * par.stretching;
					delta(i) = delta(i-1) * stretching;
				end
				XB = cumsum(delta);
			end
        
        case 'sigmoid'
        
            lambda = 1;
            
            for i = 1:par.n
                xb(i) = -10+12*(i-1)/(par.n-1);
                delta(i) = 1/(1+exp(-lambda*xb(i)));
            end
            
            delta = (delta*par.stretching+1);
            
            XB = cumsum(delta);
        
    end

end
