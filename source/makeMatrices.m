function matrices = makeMatrices(mesh,domain,boundary,numMethods) %#ok<INUSL>

% This function creates the matrices that will be used to compute derivatives and filters
% One pair of LHS and RHS will be computed for each type of region in the domain

%% Find uniform regions
findDerivativeRegions

%% Get finite differences coefficients
[centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS] = finiteDifferenceCoefficients(numMethods.spatialDerivs);
[centeredStencilLHSb, centeredStencilRHSb, decenteredStencilLHSb, decenteredStencilRHSb] = finiteDifferenceCoefficients(numMethods.spatialDerivsBuffer);
[filterStencilLHS, filterStencilRHS, filterDecenteredStencilLHS, filterDecenteredStencilRHS] = spatialFilterCoefficients(numMethods.spatialFilterStrength,numMethods.filterBorders);

%% If the Z direction is not unitary but too small to fit the filter stencil, change the stencil just for it
if mesh.nz > 1 && mesh.nz < 2*length(filterStencilRHS)-1
	filterStencilLHSz = 1;
	filterStencilRHSz = 1;
	filterDecenteredStencilLHSz = 1;
	filterDecenteredStencilRHSz = 1;
else
	filterStencilLHSz = filterStencilLHS;
	filterStencilRHSz = filterStencilRHS;
	filterDecenteredStencilLHSz = filterDecenteredStencilLHS;
	filterDecenteredStencilRHSz = filterDecenteredStencilRHS;
end

%% Make matrices
[LHSx, RHSx] = makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, derivStartsX, derivEndsX, mesh.nx, []);
[LHSy, RHSy] = makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, derivStartsY, derivEndsY, mesh.ny, []);
[LHSz, RHSz] = makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, derivStartsZ, derivEndsZ, mesh.nz, []);

[LHSxb, RHSxb] = makeMatricesEachDirection(centeredStencilLHSb, centeredStencilRHSb, decenteredStencilLHSb, decenteredStencilRHSb, derivStartsX, derivEndsX, mesh.nx, mesh.x.buffer);
[LHSyb, RHSyb] = makeMatricesEachDirection(centeredStencilLHSb, centeredStencilRHSb, decenteredStencilLHSb, decenteredStencilRHSb, derivStartsY, derivEndsY, mesh.ny, mesh.y.buffer);
[LHSzb, RHSzb] = makeMatricesEachDirection(centeredStencilLHSb, centeredStencilRHSb, decenteredStencilLHSb, decenteredStencilRHSb, derivStartsZ, derivEndsZ, mesh.nz, mesh.z.buffer);

[fLHSx, fRHSx] = makeMatricesEachDirection(filterStencilLHS, filterStencilRHS, filterDecenteredStencilLHS, filterDecenteredStencilRHS, derivStartsX, derivEndsX, mesh.nx, []);
[fLHSy, fRHSy] = makeMatricesEachDirection(filterStencilLHS, filterStencilRHS, filterDecenteredStencilLHS, filterDecenteredStencilRHS, derivStartsY, derivEndsY, mesh.ny, []);
[fLHSz, fRHSz] = makeMatricesEachDirection(filterStencilLHSz, filterStencilRHSz, filterDecenteredStencilLHSz, filterDecenteredStencilRHSz, derivStartsZ, derivEndsZ, mesh.nz, []);

%% Add buffer zones into the derivative matrices
if ~isfield(numMethods,'changeOrderX') || numMethods.changeOrderX
	[LHSx, RHSx] = addBufferToMatrix(LHSx, LHSxb, RHSx, RHSxb, mesh.nx, mesh.x.buffer);
end
if ~isfield(numMethods,'changeOrderY') || numMethods.changeOrderY
	[LHSy, RHSy] = addBufferToMatrix(LHSy, LHSyb, RHSy, RHSyb, mesh.ny, mesh.y.buffer);
end
if ~isfield(numMethods,'changeOrderZ') || numMethods.changeOrderZ
	[LHSz, RHSz] = addBufferToMatrix(LHSz, LHSzb, RHSz, RHSzb, mesh.nz, mesh.z.buffer);
end

%% Transform due to mesh stretching

LHSx = applyMetric(LHSx,mesh.X,mesh.x,domain.xf);
LHSy = applyMetric(LHSy,mesh.Y,mesh.y,domain.yf);
LHSz = applyMetric(LHSz,mesh.Z,mesh.z,domain.zf);

%% Remove ends of filters if needed
if numMethods.filterBorders
	if isfield(numMethods,'filterBordersStartX') && ~numMethods.filterBordersStartX
		for i = 1:length(fLHSx)
            fLHSx{i}(1:5,:) = 0;
            fRHSx{i}(1:5,:) = 0;
			fLHSx{i}(1:5,1:5) = eye(5);
			fRHSx{i}(1:5,1:5) = eye(5);
		end
	end
	if isfield(numMethods,'filterBordersEndX') && ~numMethods.filterBordersEndX
		for i = 1:length(fLHSx)
            fLHSx{i}(end-4:end,:) = 0;
            fRHSx{i}(end-4:end,:) = 0;
			fLHSx{i}(end-4:end,end-4:end) = eye(5);
			fRHSx{i}(end-4:end,end-4:end) = eye(5);
		end
	end
	if isfield(numMethods,'filterBordersStartY') && ~numMethods.filterBordersStartY
		for i = 1:length(fLHSy)
            fLHSy{i}(1:5,:) = 0;
            fRHSy{i}(1:5,:) = 0;
			fLHSy{i}(1:5,1:5) = eye(5);
			fRHSy{i}(1:5,1:5) = eye(5);
		end
	end
	if isfield(numMethods,'filterBordersEndY') && ~numMethods.filterBordersEndY
		for i = 1:length(fLHSy)
            fLHSy{i}(end-4:end,:) = 0;
            fRHSy{i}(end-4:end,:) = 0;
			fLHSy{i}(end-4:end,end-4:end) = eye(5);
			fRHSy{i}(end-4:end,end-4:end) = eye(5);
		end
	end
	if isfield(numMethods,'filterBordersStartZ') && ~numMethods.filterBordersStartZ
		for i = 1:length(fLHSz)
            fLHSz{i}(1:5,:) = 0;
            fRHSz{i}(1:5,:) = 0;
			fLHSz{i}(1:5,1:5) = eye(5);
			fRHSz{i}(1:5,1:5) = eye(5);
		end
	end
	if isfield(numMethods,'filterBordersEndZ') && ~numMethods.filterBordersEndZ
		for i = 1:length(fLHSz)
            fLHSz{i}(end-4:end,:) = 0;
            fRHSz{i}(end-4:end,:) = 0;
			fLHSz{i}(end-4:end,end-4:end) = eye(5);
			fRHSz{i}(end-4:end,end-4:end) = eye(5);
		end
	end
end


%% Save to output structure

matrices.x.types = typeMapX;
matrices.y.types = typeMapY;
matrices.z.types = typeMapZ;

matrices.x.LHS = LHSx;
matrices.x.RHS = RHSx;
matrices.y.LHS = LHSy;
matrices.y.RHS = RHSy;
matrices.z.LHS = LHSz;
matrices.z.RHS = RHSz;

matrices.x.fLHS = fLHSx;
matrices.x.fRHS = fRHSx;
matrices.y.fLHS = fLHSy;
matrices.y.fRHS = fRHSy;
matrices.z.fLHS = fLHSz;
matrices.z.fRHS = fRHSz;

end

function [LHS, RHS] = makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, derivStarts, derivEnds, n, bufferInfo)

    % Check if the mesh is large enough for the stencil
    if n~=1 && n < 2*length(centeredStencilRHS)-1
        error('Mesh is not large enough for one of the stencils. It has %d nodes but the stencil needs at least %d.',n,2*length(centeredStencilRHS)-1);
    end

    nTypes = length(derivStarts);
    
    LHS = cell(1,nTypes);
    RHS = cell(1,nTypes);
    
    if n == 1 % If single point in this direction, set derivative to zero and filter to one
        if centeredStencilRHS(1) == 0
            for i = 1:nTypes
                LHS{i} = 1;
                RHS{i} = 0;
            end
        else
            for i = 1:nTypes
                LHS{i} = 1;
                RHS{i} = 1;
            end
        end
        return
    end
    
    % Create the base stencils
    LHS_base = diag(centeredStencilLHS(1)*ones(n,1));
    RHS_base = diag(centeredStencilRHS(1)*ones(n,1));
    
    if centeredStencilRHS(1) == 0 % If the RHS centered stencil has center 0, its left side should have the signal inverted as this is the stencil for derivatives, otherwise, its a stencil for filters
        invertStencil = -1;
    else
        invertStencil = 1;
    end
    
    for i = 1:length(centeredStencilLHS)-1
        inds = sub2ind([n n], 1:n, 1+mod((1:n)-1+i,n));
        LHS_base(inds) = centeredStencilLHS(i+1);
        
        inds = sub2ind([n n], 1:n, 1+mod((1:n)-1-i,n));
        LHS_base(inds) = centeredStencilLHS(i+1);
    end
    
    for i = 1:length(centeredStencilRHS)-1
        inds = sub2ind([n n], 1:n, 1+mod((1:n)-1+i,n));
        RHS_base(inds) = centeredStencilRHS(i+1);
        
        inds = sub2ind([n n], 1:n, 1+mod((1:n)-1-i,n));
        RHS_base(inds) = invertStencil*centeredStencilRHS(i+1);
    end
    
    % Check if this is a buffer zone that needs to be upwind
    % and add that to the list of starts and ends so that the decentered
    % stencil is used
    if ~isempty(bufferInfo)
        if isfield(bufferInfo.i,'upwind') && bufferInfo.i.upwind
            for i = 1:nTypes
                derivStarts{i} = unique([derivStarts{i} 1:bufferInfo.i.n],'sorted');
				derivStarts{i} = derivStarts{i}(end:-1:1);
            end
        end
        if isfield(bufferInfo.f,'upwind') && bufferInfo.f.upwind
            for i = 1:nTypes
                derivEnds{i} = unique([derivEnds{i} n+1-(1:bufferInfo.f.n)],'sorted');
            end
        end
    end
    
    % Add startings and endings
    
    [mLHS, nLHS] = size(decenteredStencilLHS);
    [mRHS, nRHS] = size(decenteredStencilRHS);
    
    for i = 1:nTypes
        LHS_temp = LHS_base;
        RHS_temp = RHS_base;
        
        for j = 1:length(derivStarts{i})
            ind_start = derivStarts{i}(j);
            
            LHS_temp(ind_start:ind_start+mLHS-1,:) = 0;
            RHS_temp(ind_start:ind_start+mRHS-1,:) = 0;
            
            LHS_temp(ind_start:ind_start+mLHS-1,ind_start:ind_start+nLHS-1) = decenteredStencilLHS;
            RHS_temp(ind_start:ind_start+mRHS-1,ind_start:ind_start+nRHS-1) = decenteredStencilRHS;
            
        end
        
        for j = 1:length(derivEnds{i})
            ind_end = derivEnds{i}(j);
            
            LHS_temp(ind_end-mLHS+1:ind_end,:) = 0;
            RHS_temp(ind_end-mRHS+1:ind_end,:) = 0;
            
            LHS_temp(ind_end-mLHS+1:ind_end,ind_end-nLHS+1:ind_end) = decenteredStencilLHS(end:-1:1,end:-1:1);
            RHS_temp(ind_end-mRHS+1:ind_end,ind_end-nRHS+1:ind_end) = invertStencil*decenteredStencilRHS(end:-1:1,end:-1:1);
            
        end
        
        LHS{i} = sparse(LHS_temp);
        RHS{i} = sparse(RHS_temp);
        
    end
end

function [newMatrixL, newMatrixR] = addBufferToMatrix(baseMatrixL, bufferMatrixL, baseMatrixR, bufferMatrixR, n, bufferInfo)

ni = bufferInfo.i.n;
nf = bufferInfo.f.n;

if isfield(bufferInfo.i,'transition')
    nti = round(bufferInfo.i.transition*ni);
else
    nti = ni;
end
if isfield(bufferInfo.f,'transition')
    ntf = round(bufferInfo.f.transition*nf);
else
    ntf = nf;
end

ni1 = ni-nti+1;
ni2 = ni;
nf1 = n-nf+1;
nf2 = nf1+ntf-1;

% The buffer zone can be computed by a different type of derivatives. The transition is done smoothly.

eta = ones(n,1);

if ni > 0
    eta(1:ni1) = 0;
    eta(ni1:ni2) = 0.5-0.5*cos(linspace(0,pi,ni2-ni1+1)');
end

if nf > 0
    eta(nf2:n) = 0;
    eta(nf1:nf2) = 0.5+0.5*cos(linspace(0,pi,nf2-nf1+1)');
end

nTypes = length(baseMatrixL);
newMatrixL = cell(1,nTypes);
newMatrixR = cell(1,nTypes);

eta = sqrt(eta);

for i = 1:length(baseMatrixL)
    newMatrixL{i} = bsxfun(@times,eta,baseMatrixL{i}) + bsxfun(@times,(1-eta),bufferMatrixL{i});
    newMatrixR{i} = bsxfun(@times,eta,baseMatrixR{i}) + bsxfun(@times,(1-eta),bufferMatrixR{i});
end

end

function LHS = applyMetric(LHS,X,meshInfo,xf)

    % Scale the LHS matrix due to the mesh stretching. SL4 method is used here

    if meshInfo.periodic && meshInfo.fixPeriodicDomainSize && length(X) > 1 % Add temp node for a periodic mesh
        X = [X xf];
        addedTempNode = true;
    else
        addedTempNode = false;
    end

    n = length(X);
    if n == 1
        return
    end
    [centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS] = finiteDifferenceCoefficients('SL4');
    [LHS_temp, RHS_temp] = makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, {1}, {n}, n, []);
    
    dXdEta = LHS_temp{1}\(RHS_temp{1}*X');
    
    if addedTempNode
        dXdEta(end) = [];
    end
    
    nTypes = length(LHS);
    for i = 1:nTypes
        LHS{i} = bsxfun(@times,LHS{i},dXdEta');
    end

end
