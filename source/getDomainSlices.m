function slices = getDomainSlices(n,p)

% This function computes the size of each slice the same way the 2decomp library does
% It is used to correctly distribute the info across the processes
% The distribution is as even as possible. If an uneven distribution is needed, extra nodes are placed first in the last slices
% For example, if 10 nodes are divided across 3 slices, the division would be 3 3 4

if n == 1
    slices = [1;1];
    return
end

nPointsBase = floor(n/p);
nCeil = n - nPointsBase*p;

nPoints = ones(1,p) * nPointsBase;

nPoints(end-nCeil+1:end) = nPointsBase + 1;

slices(2,:) = cumsum(nPoints);
slices(1,1) = 1;
slices(1,2:end) = slices(2,1:end-1)+1;

end
