function p = elasticity_sampling_dx_likelihood_pars(network,delta_E,delta_S,delta_J,sigma);

% p = elasticity_sampling_dx_likelihood_pars(network,delta_E,delta_S,delta_J);
%
% parameters for elasticity_sampling_dx_likelihood.m

p.type    = 'differential expression';
p.delta_E = delta_E; 
p.delta_S = delta_S; 
p.delta_J = delta_J; 
p.network = network;
p.sigma   = sigma;