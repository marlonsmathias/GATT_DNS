function [centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS] = finiteDifferenceCoefficients(method)

% This functions contains the coefficients for each finite differences method. New methods can be defiend here if needed.

% Centered Stencils start at the center point and move to the right hand side
% Uncentered stencils are as they should be in the left hand side top of the matrix

switch method
        
	case 'SL6'

		w = 1.8;

		% Center
		matriz_a = [1	   1	      1			 -2;
					1	   4	      9			 -6;
					1	   16	      81		 -10;
					sin(w) sin(2*w)/2 sin(3*w)/3 -2*w*cos(w)];

		matriz_b = [1; 0; 0; w];

		coeffs = matriz_a\matriz_b;

		a     = coeffs(1)/2;
        b     = coeffs(2)/4;
		c	  = coeffs(3)/6;	        
		alpha = coeffs(4);

		centeredStencilLHS = [1 alpha];
        centeredStencilRHS = [0 a b c];
		
		%
		%
		% Border
		% First stage
		matriz_a = [1	1	1	1	1		1		0;
					0	1	2	3	4		5		-1;
					0	1	4	9	16		25		-2;
					0	1	8	27	64		125		-3;
					0	1	16	81	256		625		-4;
					0	1	32	243	1024	3125	-5;
					0	1	64	729	4096	15625	-6];
		
		matriz_b = [0; 1; 0; 0; 0; 0; 0]; %[0; -1; 0; 0; 0; 0; 0];

		coeffs = matriz_a\matriz_b;

		a     = coeffs(1);
        b     = coeffs(2);
		c	  = coeffs(3);
		d     = coeffs(4);
        e     = coeffs(5);
		f	  = coeffs(6);	        
		alpha = coeffs(7);
		
		matriz_lhs_aux = [1 alpha 0 0];
		matriz_rhs_aux = [a b c d e f];

		% Second stage
		matriz_a = [-1	2;
					-1	6];

		matriz_b = [-1; 0];
		
		coeffs = matriz_a\matriz_b;

		a     = coeffs(1)/2;
		alpha = coeffs(2);

		matriz_lhs_aux2 = [alpha 1 alpha 0];
		matriz_rhs_aux2 = [-a 0 a 0 0 0];

		% Third stage (same as centered SL4)

		
		matriz_a = [1      0          -2/3       ;
                    0      1          -4/3       ;
                    sin(w) sin(2*w)/2 -2*w*cos(w)];
        
        matriz_b = [4/3; -1/3; w];
        
        coeffs = matriz_a\matriz_b;
        
        a     = coeffs(1)/2;
        b     = coeffs(2)/4;
        alpha = coeffs(3);
        
        matriz_lhs_aux3 = [0 alpha 1 alpha];
        matriz_rhs_aux3 = [-b -a 0 a b 0];

		
		decenteredStencilLHS = [matriz_lhs_aux;
								matriz_lhs_aux2;
								matriz_lhs_aux3];         
        
        decenteredStencilRHS = [matriz_rhs_aux;
								matriz_rhs_aux2;
								matriz_rhs_aux3];

	case 'SL6O3'
		% Inner-Scheme: SL6 - Optimized, alpha determined minimizing the integral of dispersion
		%                     from w = [0,wf], i.e., alpha = f(integral(Edisp)|_0^wf)), wf = pi/2;
		% P3: SL4 - Optimized following the procedure above, wf = pi/2;
		% P2: C4 - Pade Scheme
		% P1: C3 -  Optimized, max. dissipation error Ediss(w=pi) = pi/2; (Edisp(w=pi/2)~6%)
		coeffs =[0.392465753424658, 1.565410958904110, 0.237260273972603, -0.017739726027397];
		alpha = coeffs(1); a = coeffs(2); b = coeffs(3); c = coeffs(4);
		centeredStencilLHS  = [1 , alpha];
		centeredStencilRHS  = [0 , a/2, b/4, c/6];

		coeffs_P3 = [0.350978473581213, 1.567318982387476, 0.134637964774951];
		alpha_P3 = coeffs_P3(1); a_P3 = coeffs_P3(2); b_P3 = coeffs_P3(3);

		decenteredStencilLHS = [    1  , (3*pi + 40)/(3*pi + 8),   0   ,     0     ;...
								  1/4 ,            1          ,  1/4   ,     0     ;...
								   0  ,        alpha_P3       ,   1   ,  alpha_P3];
		decenteredStencilRHS = [ -(13*pi + 56)/(2*(3*pi + 8)), (15*pi + 8)/(2*(3*pi + 8)), -(3*pi - 56)/(2*(3*pi + 8)), (pi - 8)/(2*(3*pi + 8)), 0;...
									 -3/4  ,    0    ,   3/4  ,     0         0   ;...
								  -b_P3/4  , -a_P3/2 ,    0   ,  a_P3/2 ,  b_P3/4 ];
										  
    case 'SL4' % 4th order spectral-like compact derivatives
        
        w = 1.8;
		           
        matriz_a = [1      0          -2/3       ;
                    0      1          -4/3       ;
                    sin(w) sin(2*w)/2 -2*w*cos(w)];
        
        matriz_b = [4/3; -1/3; w];
        
        coeffs = matriz_a\matriz_b;
        
        a     = coeffs(1)/2;
        b     = coeffs(2)/4;
        alpha = coeffs(3);
        
        centeredStencilLHS = [1 alpha];
        centeredStencilRHS = [0 a b];
        
        
        decenteredStencilLHS = [1   3 0  ;
                                1/4 1 1/4];
        
        decenteredStencilRHS = [-17/6 3/2 3/2 -1/6;
                                -3/4  0   3/4 0   ];
            
        
    case 'EX2' % 2nd order explicit finite differences
    
        centeredStencilLHS = 1;
        centeredStencilRHS = [0 1/2];
        
        decenteredStencilLHS = 1;
        decenteredStencilRHS = [-3/2 2 -1/2];
        
    case 'EX4' % 4nd order explicit finite differences
    
        centeredStencilLHS = 1;
        centeredStencilRHS = [0 2/3 -1/12];
        
        decenteredStencilLHS = 1;
        decenteredStencilRHS = [-25/12 4 -3 4/3 -1/4
                                -1/2   0 1/2 0    0 ];
        
    otherwise
        error('Finite differences method not implemented: %s. Check finiteDifferenceCoefficients file for available methods',method)
end

end
