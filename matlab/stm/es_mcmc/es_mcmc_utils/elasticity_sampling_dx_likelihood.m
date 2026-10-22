function likelihood = elasticity_sampling_dx_likelihood(result,p)

% likelihood = elasticity_sampling_dx_likelihood(result,p)
%
% Likelihood function scoring whether the response to differential
% enzyme expression (in delta_E) matches the experimental metabolite
% and flux changes (in delta_S and delta_J).
%
% result: output of es_sample_model (computed with es_options.flag_second_order = 1)
% p:      parameters, see elasticity_sampling_dx_likelihood_pars
%
% delta_E: matrix of relative enzyme changes (dE/E)
% rows: enzymes; columns: different expression patterns
% matrices delta_S (met. concentrations) and delta_J (fluxes)
% are defined accordingly
% in delta_E, all values must be given; in delta_S and delta_J,
% missing values (NaN) are allowed
%
% The predicted changes are computed by a second-order expansion
% (logarithmic), see es_simulate_perturbation

N       = p.network.N;
ind_ext = find(p.network.external);
u_ratio = 1 + p.delta_E;
s_ratio = ones(length(ind_ext),size(u_ratio,2));

r = es_simulate_perturbation(N,p.network.regulation_matrix,ind_ext,u_ratio,s_ratio,result,'logarithmic',2);

this_delta_S = r.c_ratio - 1;
this_delta_J = r.v_ratio - 1;

ind_S_finite = find(isfinite(p.delta_S));
ind_J_finite = find(isfinite(p.delta_J));

likelihood = exp( - 0.5 * (...
      sum( ( this_delta_S(ind_S_finite)-p.delta_S(ind_S_finite) ) .^2) ...
    + sum( ( this_delta_J(ind_J_finite)-p.delta_J(ind_J_finite) ) .^2) ...
        ) / p.sigma^2);
