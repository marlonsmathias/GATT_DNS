function matrix = prepareThomas(matrix)

% This function takes the matrices as inputs and prepares the vectores for the Fortran solver
% Some computations that would be repeated during runtime have their results precomputed here

% Get sizes
isPeriodic = matrix.LHS{1}(end,1) ~= 0;
nTypes = max(matrix.types(:));
N = size(matrix.LHS{1},1);

% Prepare LHS
A = zeros(N-1,nTypes);
B = zeros(N,nTypes);
C = zeros(N-1,nTypes);
D = zeros(N,nTypes);

A1 = zeros(1,nTypes);
Cn = zeros(1,nTypes);

Af = zeros(N-1,nTypes);
Bf = zeros(N,nTypes);
Cf = zeros(N-1,nTypes);
Df = zeros(N,nTypes);

A1f = zeros(1,nTypes);
Cnf = zeros(1,nTypes);

for i = 1:nTypes
    A(:,i) = full(diag(matrix.LHS{i},-1));
    B(:,i) = full(diag(matrix.LHS{i}));
    C(:,i) = full(diag(matrix.LHS{i},1));

    Af(:,i) = full(diag(matrix.fLHS{i},-1));
    Bf(:,i) = full(diag(matrix.fLHS{i}));
    Cf(:,i) = full(diag(matrix.fLHS{i},1));
end

if isPeriodic
    for i = 1:nTypes
        A1(i) = matrix.LHS{i}(1,end);
        Cn(i) = matrix.LHS{i}(end,1);

        A1f(i) = matrix.fLHS{i}(1,end);
        Cnf(i) = matrix.fLHS{i}(end,1);
    end
    
    A(1,:) = [];
    B(1,:) = [];
    C(1,:) = [];
    Af(1,:) = [];
    Bf(1,:) = [];
    Cf(1,:) = [];

    for i = 1:N-2
        C(i,:) = C(i,:)./B(i,:);
        B(i+1,:) = B(i+1,:) - A(i,:).*C(i,:);

        Cf(i,:) = Cf(i,:)./Bf(i,:);
        Bf(i+1,:) = Bf(i+1,:) - Af(i,:).*Cf(i,:);
    end

    B = 1./B;
    Bf = 1./Bf;
    
    % A1 and Cn will be stored as the first elements of A and C
    A = [A1; A];
    Af = [A1f; Af];
    C = [Cn; C];
    Cf = [Cnf; Cf];
    
    for i = 1:nTypes
        Dtemp = inv(full(matrix.LHS{i}));
        D(:,i) = Dtemp(1,:)';
        
        Dtemp = inv(full(matrix.fLHS{i}));
        Df(:,i) = Dtemp(1,:)';
    end
    
else
    A1(:) = 0;
    Cn(:) = 0;
    A1f(:) = 0;
    Cnf(:) = 0;
    
    for i = 1:N-1
        C(i,:) = C(i,:)./B(i,:);
        B(i+1,:) = B(i+1,:) - A(i,:).*C(i,:);
        
        Cf(i,:) = Cf(i,:)./Bf(i,:);
        Bf(i+1,:) = Bf(i+1,:) - Af(i,:).*Cf(i,:);
    end
    
    B = 1./B;
    Bf = 1./Bf;
    
end

% Prepare RHS

done = false;

RHSTemp = zeros(N,N,nTypes);

for i = 1:nTypes
    RHSTemp(:,:,i) = full(matrix.RHS{i});
end

RHSDiag = fullDiag(RHSTemp,0);

nDiags = 1;
while ~done
    nextDiag = fullDiag(RHSTemp,nDiags);
    prevDiag = fullDiag(RHSTemp,-nDiags);
    
    if any(nextDiag(:)) || any(prevDiag(:)) && ~(2*nDiags-1 > N)
        nDiags = nDiags + 1;
        RHSDiag = [prevDiag RHSDiag nextDiag]; %#ok<AGROW>
    else
        done = true;
    end
end

done = false;

RHSTemp = zeros(N,N,nTypes);

for i = 1:nTypes
    RHSTemp(:,:,i) = full(matrix.fRHS{i});
end

RHSDiagf = fullDiag(RHSTemp,0);

nDiagsf = 1;
while ~done
    nextDiag = fullDiag(RHSTemp,nDiagsf);
    prevDiag = fullDiag(RHSTemp,-nDiagsf);
    
    if (any(nextDiag(:)) || any(prevDiag(:))) && ~(2*nDiagsf-1 > N)
        nDiagsf = nDiagsf + 1;
        RHSDiagf = [prevDiag RHSDiagf nextDiag]; %#ok<AGROW>
    else
        done = true;
    end
end

% Store outputs
matrix.periodic = full(isPeriodic);
matrix.A = A;
matrix.B = B;
matrix.C = C;
matrix.D = D;
matrix.Af = Af;
matrix.Bf = Bf;
matrix.Cf = Cf;
matrix.Df = Df;
matrix.R = RHSDiag;
matrix.nRHS = nDiags;
matrix.Rf = RHSDiagf;
matrix.nRHSf = nDiagsf;
matrix.nTypes = nTypes;

end

function D = fullDiag(M,k)

    D = zeros(size(M,1),1,size(M,3));

    for i = 1:size(M,3)
        D(:,1,i) = diag(circshift(M(:,:,i),-k,2));
    end

end
