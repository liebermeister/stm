function E = unscale_elasticities(E, vplus, vminus, u, c);

% UNSCALE_ELASTICITIES - Convert scaled into unscaled elasticities
% 
% E = unscale_elasticities(E, vplus, vminus, u, c);

nr = length(vplus);

v  = vplus - vminus;


% -----------------------------------------------------------------------------------
% first order: simple multiply with v/p

E.un_E_plus_u   = diag(vplus)  * E.sc_E_plus_u  * diag(1./u);
E.un_E_minus_u  = diag(vminus) * E.sc_E_minus_u * diag(1./u)  ;
E.un_E_plus_c   = diag(vplus)  * E.sc_E_plus_c  * diag(1./c)  ;
E.un_E_minus_c  = diag(vminus) * E.sc_E_minus_c * diag(1./c)  ;
	        
E.un_E_u        = E.un_E_plus_u - E.un_E_minus_u;
E.un_E_c        = E.un_E_plus_c - E.un_E_minus_c;

E.sc_E_u        = diag(1./v) * E.un_E_u * diag(u);
E.sc_E_c        = diag(1./v) * E.un_E_c * diag(c);


% -----------------------------------------------------------------------------------
% second order 

E.un_E_plus_cc   = tensor_scale(E.sc_E_plus_cc,1,vplus);
E.un_E_plus_cc   = tensor_scale(E.un_E_plus_cc,2,1./c);
E.un_E_plus_cc   = tensor_scale(E.un_E_plus_cc,3,1./c);
	        
E.un_E_minus_cc  = tensor_scale(E.sc_E_minus_cc,1,vminus);
E.un_E_minus_cc  = tensor_scale(E.un_E_minus_cc,2,1./c);
E.un_E_minus_cc  = tensor_scale(E.un_E_minus_cc,3,1./c);
	        
E.un_E_cc        = E.un_E_plus_cc - E.un_E_minus_cc;

E.sc_E_cc        = tensor_scale(E.un_E_cc,1,1./v);
E.sc_E_cc        = tensor_scale(E.sc_E_cc,2,c);
E.sc_E_cc        = tensor_scale(E.sc_E_cc,3,c);


% enzyme elasiticities

E.un_E_plus_cu   = tensor_scale(E.sc_E_plus_cu,1,vplus);
E.un_E_plus_cu   = tensor_scale(E.un_E_plus_cu,2,1./c);
E.un_E_plus_cu   = tensor_scale(E.un_E_plus_cu,3,1./u);

E.un_E_minus_cu  = tensor_scale(E.sc_E_minus_cu,1,vminus);
E.un_E_minus_cu  = tensor_scale(E.un_E_minus_cu,2,1./c);
E.un_E_minus_cu  = tensor_scale(E.un_E_minus_cu,3,1./u);

E.un_E_cu        = E.un_E_plus_cu - E.un_E_minus_cu;

E.sc_E_cu        = tensor_scale(E.un_E_cu,1,1./v);
E.sc_E_cu        = tensor_scale(E.sc_E_cu,2,c);
E.sc_E_cu        = tensor_scale(E.sc_E_cu,3,u);

E.un_E_plus_uu   = zeros(nr,nr,nr);
E.un_E_minus_uu  = zeros(nr,nr,nr);
E.un_E_uu        = zeros(nr,nr,nr);
E.sc_E_plus_uu   = zeros(nr,nr,nr);
E.sc_E_minus_uu  = zeros(nr,nr,nr);
E.sc_E_uu        = zeros(nr,nr,nr);
