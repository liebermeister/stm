function p = elasticity_sampling_tc_likelihood_pars(network,c_init,data);

% p = elasticity_sampling_tc_likelihood_pars(network,c_init,data);
%
% parameters for elasticity_sampling_tc_likelihood 

p.type    = 'timecourse';
p.data    = data; 
p.c_init  = c_init; 
p.network = network;
