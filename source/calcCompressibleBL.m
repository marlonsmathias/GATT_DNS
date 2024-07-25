function flow = calcCompressibleBL(flowParameters,adiabWall,mesh)
    Re      = flowParameters.Re;       % Reynolds Number
    Minf    = flowParameters.Ma;       % Mach Number
    Pr      = flowParameters.Pr;       % Prandtl Number
    Tinf    = flowParameters.T0;       % Temp. [K]
    gamma   = flowParameters.gamma;    % Gamma

    
    xR      = Re / 1.7208^2; % Reference x (where delta*(xR) = 1, assuming imcompressible Blasius)
    E0      = 1/(gamma*(gamma-1)*Minf^2);
    
    Twall   = 1; % Twall/Tinf
    if isfield(flowParameters,'Twall')
        Twall = flowParameters.Twall / Tinf; % Twall/Tinf
    end
    sol = solve_compressibleBL(Minf,Pr,Tinf,Twall,gamma,adiabWall);

    eta_xR = sol.eta;
    U_xR   = sol.f_p;
    R_xR   = 1 ./ sol.rbar;
    E_xR   = sol.rbar * E0;
    V_xR   = -(sol.rbar).*(sol.f - eta_xR' .* sol.f_p)*(1.7208/sqrt(2))*(1/Re);
    
    y_xR = eta_xR*sqrt(2)/1.7208;
    dS_xR = trapz(y_xR',1-U_xR .* R_xR / (U_xR(end) * R_xR(end)));
    
    fprintf(' Initial Flow: Compressible BL\n');
    if adiabWall
        fprintf('   Adiabatic Wall, (dT/dy)_wall = 0, (Twall/TInf) = %1.3f \n',E_xR(1)/E0);
    else        
        fprintf('   Isothermal Wall, T_wall = cte,    (Twall/TInf) = %1.3f \n',E_xR(1)/E0);
    end

    fprintf('   !!! BL thickness, (effective) / (incomp. BL, Blasius) = %1.4f\n',dS_xR);


    %% --------------------------
    X = mesh.X; nx = length(X);
    Y = mesh.Y; ny = length(Y);
    Z = mesh.Z; nz = length(Z);

    R = ones(nx,ny);
    U = ones(nx,ny);
    V = zeros(nx,ny);
    W = zeros(nx,ny);
    E = ones(nx,ny) * E0;

       
    indX = find(X>0);
    indY = find(Y>=0);
    for ix = 1:length(indX)        
        y_xL = y_xR*sqrt(X(indX(ix)) / xR);

        U(indX(ix),indY) = interp1(y_xL,U_xR',Y(indY),'makima',U_xR(end));
        V(indX(ix),indY) = interp1(y_xL,V_xR',Y(indY),'makima',V_xR(end));
        R(indX(ix),indY) = interp1(y_xL,R_xR',Y(indY),'makima',R_xR(end));
        E(indX(ix),indY) = interp1(y_xL,E_xR',Y(indY),'makima',E_xR(end));
        clearvars y_xL
    end
    clearvars indY

    U(X==0,Y==0) = 0;
    indY = find(Y<0);
    U(:,indY) = repmat(U(:,Y==0),[1 length(indY)]);
    V(:,indY) = repmat(V(:,Y==0),[1 length(indY)]);
    R(:,indY) = repmat(R(:,Y==0),[1 length(indY)]);
    E(:,indY) = repmat(E(:,Y==0),[1 length(indY)]);
    
    U = repmat(U,[1 1 nz]);
    V = repmat(V,[1 1 nz]);
    W = repmat(W,[1 1 nz]);
    R = repmat(R,[1 1 nz]);
    E = repmat(E,[1 1 nz]);

    flow.U = U;
    flow.V = V;
    flow.W = W;
    flow.R = R;
    flow.E = E;
end

function sol = solve_compressibleBL(Minf,Pr,Tinf,Twall,Gamma,adiabWall)
    
    C2  = 110;     % Sutherland Coefficient [Kelvin]
    lim = 10;       % The value which simulates lim-> inf
    N   = 500;      % Number of Point
    h   = lim/N;      % Delta y 
    delta = 1e-10;   % Small Number for shooting method
    eps = 1e-9;
    
    adi = 0;
    if adiabWall
        adi = 1;
    end

    % Initializing
    y1    = zeros(N+1,1);   % f
    y2    = zeros(N+1,1);   % f'
    y3    = zeros(N+1,1);   % f''
    y4    = zeros(N+1,1);   % rho(eta)
    y5    = zeros(N+1,1);   % rho(eta)'
    eta   = 0:h:lim;        % Iteration of eta up to infinity
    dalfa = 0;
    dbeta = 0;

    if adi==1
        % Boundary Conditions for Adiabatic Case
        y1(1) = 0;
        y2(1) = 0;
        y5(1) = 0;

        % Initial Guess for the beginning of simulation
        alfa0 = 0.1;      % Initial Guess
        beta0 = 3;        % Initial Guess
    elseif adi==0
        % Boundary Conditions for Isothermal Case
        y1(1) = 0;
        y2(1) = 0;
        y4(1) = Twall;

        % Initial Guess for Beginning of Simulation
        alfa0 = 0.1;      % Initial Guess
        beta0 = 3;        % Initial Guess  
    end
    
    for ite = 1:100000 
        if adi==1
            % Boundary Conditions for Adiabatic Case
            y1(1) = 0;
            y2(1) = 0;
            y5(1) = 0;

            y3(1) = alfa0;
            y4(1) = beta0;
        elseif adi==0
            % Boundary Conditions for Isothermal Case
            y1(1) = 0;
            y2(1) = 0;
            y4(1) = Twall;

            y3(1) = alfa0;
            y5(1) = beta0;
        end
        [y1,y2,y3,y4,y5] = RK(eta,h,y1,y2,y3,y4,y5,C2,Tinf,Minf,Pr,Gamma);

        y2old = y2(end);
        y4old = y4(end);

        if adi==1
            % Boundary Conditions for Adiabatic Case
            y1(1) = 0;
            y2(1) = 0;
            y5(1) = 0;

            y3(1) = alfa0+delta;
            y4(1) = beta0;
        elseif adi==0
            % Boundary Conditions for Isothermal Case
            y1(1) = 0;
            y2(1) = 0;
            y4(1) = Twall;

            y3(1) = alfa0+delta;
            y5(1) = beta0;
        end
        [y1,y2,y3,y4,y5] = RK(eta,h,y1,y2,y3,y4,y5,C2,Tinf,Minf,Pr,Gamma);

        y2new1 = y2(end);
        y4new1 = y4(end);

        if adi==1
            % Boundary Conditions for Adiabatic Case
            y1(1) = 0;
            y2(1) = 0;
            y5(1) = 0;

            y3(1) = alfa0;
            y4(1) = beta0+delta;
        elseif adi==0
            % Boundary Conditions for Isothermal Case
            y1(1) = 0;
            y2(1) = 0;
            y4(1) = Twall;

            y3(1) = alfa0;
            y5(1) = beta0+delta;
        end
        [y1,y2,y3,y4,y5] = RK(eta,h,y1,y2,y3,y4,y5,C2,Tinf,Minf,Pr,Gamma);

        y2new2 = y2(end);
        y4new2 = y4(end);

        a11 = (y2new1-y2old)/delta;
        a21 = (y4new1-y4old)/delta;
        a12 = (y2new2-y2old)/delta;
        a22 = (y4new2-y4old)/delta;
        r1 = 1-y2old;
        r2 = 1-y4old;
        dalfa = (a22*r1-a12*r2)/(a11*a22-a12*a21);
        dbeta = (a11*r2-a21*r1)/(a11*a22-a12*a21);
        alfa0 = alfa0 + dalfa;
        beta0 = beta0 + dbeta;

        if (abs(y2(end)-1)<eps) && (abs(y4(end)-1)<eps)
            Truey2 = y2(1);
            Truey4 = y4(1);
            break
        end
    end
    
    xaxis = zeros(length(eta),1);
    for i=2:length(eta)
        xaxis(i) = (eta(i)-0)*(y4(1)+2*sum(y4(2:i-1))+y4(i))/(2*eta(i))*h;
    end

    sol.eta     = xaxis';
%     sol.eta     = eta;
    sol.f       = y1;
    sol.f_p     = y2;
    sol.f_pp    = y3;
    sol.rbar    = y4;
    sol.rbar_p  = y5;
end


function [y2] = Y1(y2)
end

function [y3] = Y2(y3)
end

function [RHS] = Y3(y1,y3,y4,y5,C2,Tinf)
    RHS = -y3*((y5/(2*(y4)))-(y5/(y4+C2/Tinf))) ...
                         -y1*y3*((y4+C2/Tinf)/(sqrt(y4)*(1+C2/Tinf)));
end

function [y5] = Y4(y5)
end

function [RHS] = Y5(y1,y3,y4,y5,C2,Tinf,Minf,Pr,Gamma)
    RHS = -y5^2*((0.5/y4)-(1/(y4+C2/Tinf)))...
                         -Pr*y1*y5/sqrt(y4)*(y4+C2/Tinf)/(1+C2/Tinf)...
                         -(Gamma-1)*Pr*Minf^2*y3^2;
end

function [y1,y2,y3,y4,y5] = RK(eta,h,y1,y2,y3,y4,y5,C2,Tinf,Minf,Pr,Gamma)
    for i=1:(length(eta)-1)
        
        k11 = Y1(y2(i));
        k21 = Y2(y3(i));
        k31 = Y3(y1(i), y3(i), y4(i), y5(i),C2,Tinf);
        k41 = Y4(y5(i));
        k51 = Y5(y1(i), y3(i), y4(i), y5(i),C2,Tinf,Minf,Pr,Gamma);
        
        
        k12 = Y1(y2(i)+0.5*h*k21);
        k22 = Y2(y3(i)+0.5*h*k31);
        k32 = Y3(y1(i)+0.5*h*k11, y3(i)+0.5*h*k31, y4(i)+0.5*h*k41, y5(i)+0.5*h*k51,C2,Tinf);
        k42 = Y4(y5(i)+0.5*h*k51);
        k52 = Y5(y1(i)+0.5*h*k11, y3(i)+0.5*h*k31, y4(i)+0.5*h*k41, y5(i)+0.5*h*k51,C2,Tinf,Minf,Pr,Gamma);
        
        k13 = Y1(y2(i)+0.5*h*k22);
        k23 = Y2(y3(i)+0.5*h*k32);
        k33 = Y3(y1(i)+0.5*h*k12, y3(i)+0.5*h*k32, y4(i)+0.5*h*k42, y5(i)+0.5*h*k52,C2,Tinf);
        k43 = Y4(y5(i)+0.5*h*k52);
        k53 = Y5(y1(i)+0.5*h*k12, y3(i)+0.5*h*k32, y4(i)+0.5*h*k42, y5(i)+0.5*h*k52,C2,Tinf,Minf,Pr,Gamma);
        
        k14 = Y1(y2(i)+h*k23);
        k24 = Y2(y3(i)+h*k33);
        k34 = Y3(y1(i)+h*k13, y3(i)+h*k33, y4(i)+h*k43, y5(i)+h*k53,C2,Tinf);
        k44 = Y4(y5(i)+h*k53);
        k54 = Y5(y1(i)+h*k13, y3(i)+h*k33, y4(i)+h*k43, y5(i)+h*k53,C2,Tinf,Minf,Pr,Gamma);
        
        y5(i+1) = y5(i) + (1/6)*(k51 + 2*k52 + 2*k53 + k54)*h;
        y4(i+1) = y4(i) + (1/6)*(k41 + 2*k42 + 2*k43 + k44)*h;
        y3(i+1) = y3(i) + (1/6)*(k31 + 2*k32 + 2*k33 + k34)*h; 
        y2(i+1) = y2(i) + (1/6)*(k21 + 2*k22 + 2*k23 + k24)*h;
        y1(i+1) = y1(i) + (1/6)*(k11 + 2*k12 + 2*k13 + k14)*h;
    end
end