Tin = [11700:16600];
pastaIn = 'rugosidade6_5_BC2';

Tout = Tin(end) + 1;
pastaOut = pastaIn;

for t = Tin
    
    tstr = sprintf('%s/flow_%010d.mat',pastaIn,t);
    t
    current = load(tstr);

	if t == Tin(1)
		U = current.U;
		V = current.V;
		W = current.W;
		R = current.R;
		E = current.E;
	else
		U = U + current.U;
		V = V + current.V;
		W = W + current.W;
		R = R + current.R;
		E = E + current.E;
	end
end

U = U/length(Tin);
V = V/length(Tin);
W = W/length(Tin);
R = R/length(Tin);
E = E/length(Tin);
t = current.t;

tstr = sprintf('%s/flow_%010d.mat',pastaOut,Tout);
save(tstr,'t','U','V','W','R','E')
