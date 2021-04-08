function [filterStencilLHS, filterStencilRHS, filterDecenteredStencilLHS, filterDecenteredStencilRHS] = spatialFilterCoefficients(alpha,filterBorders)

% This function outputs the coefficients for the 10th order spatial filter from Gaitonde 1998 for a given alpha

filterStencilLHS = [1 alpha];

filterStencilRHS = [(193+126*alpha)/256 (105+302*alpha)/512 (-15+30*alpha)/128 (45-90*alpha)/1024 (-5+10*alpha)/512 (1-2*alpha)/1024];

if islogical(filterBorders)
	if filterBorders
		filterBorders = 'decentered';
	else
		filterBorders = 'off';
	end
end

switch filterBorders
	case 'decentered'
		filterDecenteredStencilLHS = [1 alpha];
		
		a_bound( 1,5) = (-1+2*alpha)/1024;
		a_bound( 2,5) = (5-10*alpha)/512;
		a_bound( 3,5) = (-45+90*alpha)/1024;
		a_bound( 4,5) = (15+98*alpha)/128;
		a_bound( 5,5) = (407+210*alpha)/512;
		a_bound( 6,5) = (63+130*alpha)/256;
		a_bound( 7,5) = (-105+210*alpha)/512;
		a_bound( 8,5) = (15-30*alpha)/128;
		a_bound( 9,5) = (-45+90*alpha)/1024;
		a_bound(10,5) = (5-10*alpha)/512;
		a_bound(11,5) = (-1+2*alpha)/1024;

		a_bound( 1,4) = (1-2*alpha)/1024;
		a_bound( 2,4) = (-5+10*alpha)/512;
		a_bound( 3,4) = (45+934*alpha)/1024;
		a_bound( 4,4) = (113+30*alpha)/128;
		a_bound( 5,4) = (105+302*alpha)/512;
		a_bound( 6,4) = (-63+126*alpha)/256;
		a_bound( 7,4) = (105-210*alpha)/512;
		a_bound( 8,4) = (-15+30*alpha)/128;
		a_bound( 9,4) = (45-90*alpha)/1024;
		a_bound(10,4) = (-5+10*alpha)/512;
		a_bound(11,4) = (1-2*alpha)/1024;

		a_bound( 1,3) = (-1+2*alpha)/1024;
		a_bound( 2,3) = (5+502*alpha)/512;
		a_bound( 3,3) = (979+90*alpha)/1024;
		a_bound( 4,3) = (15+98*alpha)/128;
		a_bound( 5,3) = (-105+210*alpha)/512;
		a_bound( 6,3) = (63-126*alpha)/256;
		a_bound( 7,3) = (-105+210*alpha)/512;
		a_bound( 8,3) = (15-30*alpha)/128;
		a_bound( 9,3) = (-45+90*alpha)/1024;
		a_bound(10,3) = (5-10*alpha)/512;
		a_bound(11,3) = (-1+2*alpha)/1024;

		a_bound( 1,2) = (1+1022*alpha)/1024;
		a_bound( 2,2) = (507+10*alpha)/512;
		a_bound( 3,2) = (45+934*alpha)/1024;
		a_bound( 4,2) = (-15+30*alpha)/128;
		a_bound( 5,2) = (105-210*alpha)/512;
		a_bound( 6,2) = (-63+126*alpha)/256;
		a_bound( 7,2) = (105-210*alpha)/512;
		a_bound( 8,2) = (-15+30*alpha)/128;
		a_bound( 9,2) = (45-90*alpha)/1024;
		a_bound(10,2) = (-5+10*alpha)/512;
		a_bound(11,2) = (1-2*alpha)/1024;

		a_bound( 1,1) = (1023+1*alpha)/1024;
		a_bound( 2,1) = (5+507*alpha)/512;
		a_bound( 3,1) = (-45+45*alpha)/1024;
		a_bound( 4,1) = (15-15*alpha)/128;
		a_bound( 5,1) = (-105+105*alpha)/512;
		a_bound( 6,1) = (63-63*alpha)/256;
		a_bound( 7,1) = (-105+105*alpha)/512;
		a_bound( 8,1) = (15-15*alpha)/128;
		a_bound( 9,1) = (-45+45*alpha)/1024;
		a_bound(10,1) = (5-5*alpha)/512;
		a_bound(11,1) = (-1+1*alpha)/1024;
		
		filterDecenteredStencilRHS = a_bound';
		
	case 'reducedOrder'
        filterDecenteredStencilLHS = 1;
        
		F2 = [1/2+alpha, 1/2+alpha]./[1 2];
		F4 = [5/8+3/4*alpha, 1/2+alpha, -1/8+1/4*alpha]./[1 2 2];
		F6 = [11/16+5/8*alpha, 15/32+17/16*alpha, -3/16+3/8*alpha, 1/32-1/16*alpha]./[1 2 2 2];
		F8 = [93/128+70/128*alpha, 7/16+18/16*alpha, -7/32+14/32*alpha, 1/16-1/8*alpha, -1/128+1/64*alpha]./[1 2 2 2 2];
		
		filterDecenteredStencilRHS = zeros(5,9);
		
		filterDecenteredStencilRHS(1,1) = 1;
		filterDecenteredStencilRHS(2,1:3) = F2([2 1 2]);
		filterDecenteredStencilRHS(3,1:5) = F4([3 2 1 2 3]);
		filterDecenteredStencilRHS(4,1:7) = F6([4 3 2 1 2 3 4]);
		filterDecenteredStencilRHS(5,1:9) = F8([5 4 3 2 1 2 3 4 5]);
		
	case 'off'
		filterDecenteredStencilLHS = diag(ones(1,5));
		filterDecenteredStencilRHS = diag(ones(1,5));
end

end
