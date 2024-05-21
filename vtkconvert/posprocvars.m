function [] = posprocvars(fstart,fend,p_col,p_row)
    for Flow= fstart:fend
        m = load('mesh.mat');
        
        flow = load(sprintf('flow_%010d.mat',Flow));
        flowname = strcat(sprintf('flow_%010d.mat',Flow));
        U = flow.U;
        U(isnan(U)) = 0;
        V = flow.V;
        V(isnan(V)) = 0;
        W = flow.W;
        W(isnan(W)) = 0;

        [Uy, Ux, Uz] = gradient(U);
        [Vy, Vx, Vz] = gradient(V);
        [Wy, Wx, Wz] = gradient(W);

        Ux = Ux./gradient(m.X');
        Uy = Uy./gradient(m.Y);
        Uz = Uz./permute(gradient(m.Z),[1 3 2]);

        Vx = Vx./gradient(m.X');
        Vy = Vy./gradient(m.Y);
        Vz = Vz./permute(gradient(m.Z),[1 3 2]);

        Wx = Wx./gradient(m.X');
        Wy = Wy./gradient(m.Y);
        Wz = Wz./permute(gradient(m.Z),[1 3 2]);



        L2 = 0*U;

        for i = 1:size(U,1)
            for j = 1:size(U,2)-60
                if any(U(i,j,:)~=0)
                    for k = 1:size(U,3)
                        L2(i,j,k) = lambda2([Ux(i,j,k) Uy(i,j,k) Uz(i,j,k); Vx(i,j,k) Vy(i,j,k) Vz(i,j,k); Wx(i,j,k) Wy(i,j,k) Wz(i,j,k)]);
                        if L2(i,j,k) >= 0
                            L2(i,j,k) = NaN;
                        end
                    end
                end
            end
        end

        q = 0*U;
        qs = 0*U;

        S11 = Ux;
        S12 = 0.5*(Uy+Vx);
        S13 = 0.5*(Uz+Wx);
        S22 = Vy;
        S23 = 0.5*(Vz+Wy);
        S33 = Wz;
        Omga12 = 0.5*(Uy-Vx);
        Omga13 = 0.5*(Uz-Wx);
        Omga23 = 0.5*(Vz-Wy);

        for i = 60:size(U,1)
            for j = 1:size(U,2)-60
                if any(U(i,j,:)~=0)
                    for k = 1:size(U,3)
                        q(i,j,k) = qcrit(Omga12(i,j,k),Omga13(i,j,k),Omga23(i,j,k),S11(i,j,k),S22(i,j,k),S33(i,j,k),S12(i,j,k),S13(i,j,k),S23(i,j,k));
                        qs(i,j,k) = NaN;
                        if q(i,j,k) <= 0
                            qs(i,j,k) = q(i,j,k);
                            q(i,j,k) = NaN;
                        end
                    end
                end
            end
        end
        mq = max(max(max(q)));
        Q = q./mq;
            
        Dilatation = (Ux+Vy+Wz);

        aux3 = load('mesh.mat');
        gama = 1.4;
        X = aux3.X;
        Y = aux3.Y;
        Z = aux3.Z;
        wall = aux3.wall;
        flowParameters = aux3.flowParameters;
        flowType = aux3.flowType;
        % P = 0.4.*flow.R.*flow.E;
        % T = flow.E.*gama.*(gama-1).*(flowParameters.Ma.*flowParameters.Ma);
        [Vorti,Vortj,Vortk,Vort] = vorticity(flow);
        Mach = (1/sqrt(gama.*(gama-1))).*flow.U./sqrt(flow.E);

        save(flowname,'L2','Q','Dilatation','Vortk','Mach','-append')

        clear flow aux1 aux2 aux3 L2 Q X dilatation U V W Ux Uy Uz Vx Vy Vz Wx Wy Wz mq q qs S11 S12 S13 S22 S23 S33 Omega12 Omega13 Omega23;
    
    end    
end



function [Vorti,Vortj,Vortk,Vort] = vorticity(flow)

    [Vorti,Vortj,Vortk] = curl(flow.U,flow.V,flow.W);
    Vort = sqrt(Vorti.*Vorti + Vortj.*Vortj + Vortk.*Vortk);

end

function l2 = lambda2(J)
    S = (J + J')/2;
    O = (J - J')/2;

    L = eig(S^2 + O^2);
    l2 = real(L(2));
end

function q = qcrit(Omga12,Omga13,Omga23,S11,S22,S33,S12,S13,S23)
    A = 2*(Omga12^2) + 2*(Omga13^2) + 2*(Omga23^2) - S11^2 - S22^2 - S33^2 - 2*(S12^2) - 2*(S13^2) - 2*(S23^2);
    q = 0.5*A;
end