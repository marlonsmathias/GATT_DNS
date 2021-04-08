% This script is called if the SFD is to be applied only at the buffer zones.
% It generates an array the same size as the domain which smoothly goes from 0 at the physical region to SFD_X at the buffer zones.

SFD_X = ones(mesh.nx,mesh.ny,mesh.nz);

% Set default value for applyX,Y,Z
if ~isfield(numMethods.SFD,'applyX') || (islogical(numMethods.SFD.applyX) && numMethods.SFD.applyX)
	numMethods.SFD.applyX = 2;
end
if ~isfield(numMethods.SFD,'applyY') || (islogical(numMethods.SFD.applyY) && numMethods.SFD.applyY)
	numMethods.SFD.applyY = 2;
end
if ~isfield(numMethods.SFD,'applyZ') || (islogical(numMethods.SFD.applyZ) && numMethods.SFD.applyZ)
	numMethods.SFD.applyZ = 2;
end

if numMethods.SFD.applyX == 2 || numMethods.SFD.applyX == -1
	for i = 1:mesh.x.buffer.i.n
		SFD_X(i,:,:) = SFD_X(i,:,:) * (0.5 - 0.5*cos(pi*(i-1)/(mesh.x.buffer.i.n)));
	end
end
if numMethods.SFD.applyX == 2 || numMethods.SFD.applyX == 1
	for i = 1:mesh.x.buffer.f.n
		SFD_X(end-i+1,:,:) = SFD_X(end-i+1,:,:) * (0.5 - 0.5*cos(pi*(i-1)/(mesh.x.buffer.f.n)));
	end
end

if numMethods.SFD.applyY == 2 || numMethods.SFD.applyY == -1
	for j = 1:mesh.y.buffer.i.n
		SFD_X(:,j,:) = SFD_X(:,j,:) * (0.5 - 0.5*cos(pi*(j-1)/(mesh.y.buffer.i.n)));
	end
end
if numMethods.SFD.applyY == 2 || numMethods.SFD.applyY == 1
	for j = 1:mesh.y.buffer.f.n
		SFD_X(:,end-j+1,:) = SFD_X(:,end-j+1,:) * (0.5 - 0.5*cos(pi*(j-1)/(mesh.y.buffer.f.n)));
	end
end

if numMethods.SFD.applyZ == 2 || numMethods.SFD.applyZ == -1
	for k = 1:mesh.z.buffer.i.n
		SFD_X(:,:,k) = SFD_X(:,:,k) * (0.5 - 0.5*cos(pi*(k-1)/(mesh.z.buffer.i.n)));
	end
end
if numMethods.SFD.applyZ == 2 || numMethods.SFD.applyZ == 1
	for k = 1:mesh.z.buffer.f.n
		SFD_X(:,:,end-k+1) = SFD_X(:,:,end-k+1) * (0.5 - 0.5*cos(pi*(k-1)/(mesh.z.buffer.f.n)));
	end
end

SFD_X = numMethods.SFD.X * (1-SFD_X);

if isfield(numMethods.SFD,'extraRegion')
    for i = 1:length(numMethods.SFD.extraRegion)
        ER = numMethods.SFD.extraRegion{i};
        % Get radius from center point
        R = sqrt((mesh.X'-ER.location(1)).^2/ER.size(1)^2 + (mesh.Y-ER.location(2)).^2/ER.size(2)^2 + (permute(mesh.Z,[1 3 2])-ER.location(3)).^2/ER.size(3)^2);
        R(R>1) = 1;
        R = 0.5+0.5*cos(pi*R);
        SFD_X = SFD_X + ER.X * R;
    end
    clear ER R
end
