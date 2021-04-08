function blocks = getMatrixTypeBlocks(typeMap,p_row,p_col)

% This function transforms an array of derivative types into a list of coordinates for each type
% For example, for derivatives in X, the input matrix will be ny x nz and each entry is the type of derivative to be applied at each row
% The output is a list with 5 columns. The first is the derivative type, columns 2 to 5 are the limits for each region in Y and Z
% One list is created for each MPI process

maxBlockSize = 128; % Maximum number of rows per block, larger blocks will be divided

blocks = cell(p_row*p_col,1);

[J,K] = size(typeMap);

% Get all blocks

allBlocks = [];
for k = 1:K
    starts = [1;find(diff(typeMap(:,k)) ~= 0)+1];
    ends = [starts(2:end)-1;J];
    allBlocks = [allBlocks; typeMap(starts,k) starts ends k*ones(length(starts),2)]; %#ok<*AGROW>
end

% Merge blocks
i = 1;
while i < size(allBlocks,1)
    
    j = i+1;
    while j <= size(allBlocks,1)
        
        if all(allBlocks(i,1:3) == allBlocks(j,1:3)) && allBlocks(i,5)+1 == allBlocks(j,4)
            allBlocks(i,5) = allBlocks(j,5);
            allBlocks(j,:) = [];
        else
            j = j + 1;
        end
    end
    
    i = i+1;
end

% Divide blocks for processors

domainSlicesY = getDomainSlices(J,p_row);
domainSlicesZ = getDomainSlices(K,p_col);

for j = 1:p_row
    for k = 1:p_col
        nProc = k + (j-1)*p_col;
        
        bL = allBlocks;
        
        Ji = domainSlicesY(1,j);
        Jf = domainSlicesY(2,j);
        Ki = domainSlicesZ(1,k);
        Kf = domainSlicesZ(2,k);
        
        bL(bL(:,2) > Jf | bL(:,3) < Ji | bL(:,4) > Kf | bL(:,5) < Ki,:) = [];
        
        bL(:,2) = max(bL(:,2),Ji);
        bL(:,3) = min(bL(:,3),Jf);
        bL(:,4) = max(bL(:,4),Ki);
        bL(:,5) = min(bL(:,5),Kf);
        
        blocks{nProc} = bL;
        
    end
end

% Reduce block sizes if needed
if ~isinf(maxBlockSize)
    for nProc = 1:p_col*p_row
        bL = blocks{nProc};
        bLnew = zeros(0,5);
        for j = 1:size(bL,1)
            iSize = bL(j,3)-bL(j,2)+1;
            jSize = bL(j,5)-bL(j,4)+1;
            bSize = iSize * jSize;
            nSlices = ceil(bSize/maxBlockSize);
            if nSlices == 1
                bLnew(end+1,:) = bL(j,:);
            else
                iSlices = ceil(nSlices/jSize);
                jSlices = floor(nSlices/iSlices);

                iSlicesInd = getDomainSlices(iSize,iSlices);
                jSlicesInd = getDomainSlices(jSize,jSlices);

                iSlicesInd = bL(j,2) + iSlicesInd' - 1;
                jSlicesInd = bL(j,4) + jSlicesInd' - 1;

                for ii = 1:iSlices
                    for jj = 1:jSlices
                        bLnew(end+1,:) = [bL(j,1) iSlicesInd(ii,:) jSlicesInd(jj,:)];
                    end
                end
            end
        end
        blocks{nProc} = bLnew;
    end
end

end
