function result = sample_elasticities(S,R,ind_ext,v,c,keq,options)
  
% result = sample_elasticities()
%
% Convenience (wrapper) function for es_sample_elasticities with default values
%
% Elasticity sampling with thermodynamic forces 
% The sampling of elasticities is based on a sampling of saturation values (called beta values)
%
% Input:
%   S        Stoichiometric matrix (size: metabolites x reactions)
%            note that the matrix contains ALL metabolites (not just internal metabolites)
%   R        Regulation matrix     (size: reactions x metabolites)
%            (elements: 1 for activation, -1 for inhibition, otherwise 0)
%   ind_ext  index vector of external metabolites (currently not used by the code)
%   v        Fluxes (vector)
%   c        Metabolite concentrations (vector)
%   keq      Equilibrium constants (vector)
%            Note that keq can be easily computed from a vector dg0 of Delta G0 values: keq = exp(-dg0/RT)
%            or from a vector of mu0 values: keq = exp(-1/RT * N' * mu0)
%   options  Currently not used
  
[es_options, es_constraints] = es_default_options(S);
c0 = exp(pinv(log(keq)));
u = ones(size(S,2),1);
h = ones(size(S,2));
A = RT * [log(keq)-S'*log(c)];
result = es_sample_elasticities(S, R, ind_ext, es_constraints, es_options, c0, c, u, v, keq, [],[],[],[],A);
