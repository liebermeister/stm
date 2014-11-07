function likelihood = elasticity_sampling_dx_likelihood(result,p,kinetic_law)

% likelihood = elasticity_sampling_dx_likelihood(result,p, kinetic_law)
%
% Likelihood function scoring whether the response to differential
% enzyme expression (in delta_E) matches the experimental metabolite 
% and flux changes (in delta_S and delta_J).
%
% delta_E: matrix of relative enzyme changes (dE/E)
% rows: enzymes; columns: different expression patterns
% matrices delta_S (met. concentrations) and delta_J (fluxes)
% are defined accordingly
% in delta_E, all values must be given; in delta_S and delta_J,
% missing values (NaN) are allowed

% EXAMPLE: kinetic_law = 'cs'; 

evela(default('kinetic_law','''cs'''));

p.network.kinetics =  set_kinetics(p.network,kinetic_law,result.parameters);
 
[this_delta_S,this_delta_J] = predict_response_to_enzyme_change(...
     p.delta_E,result.control.RSu_sc,result.control.RJu_sc,result.control.RSuu_sc,result.control.RJuu_sc);

ind_S_finite = find(isfinite(p.delta_S));
ind_J_finite = find(isfinite(p.delta_J));

likelihood = exp( - 0.5 * (...
      sum( ( this_delta_S(ind_S_finite)-p.delta_S(ind_S_finite) ) .^2) ...
    + sum( ( this_delta_J(ind_J_finite)-p.delta_J(ind_J_finite) ) .^2) ...
        ) / p.sigma^2);
