function likelihood = elasticity_sampling_likelihood(result,p)

% likelihood = elasticity_sampling_likelihood(result,p)

parameters.u       = result.u;     
parameters.c       = result.c;     
parameters.KA      = result.KA;    
parameters.KI      = result.KI;    
parameters.KM      = result.KM;    
parameters.KV      = result.KV;    
parameters.Keq = result.Keq;

p.network.kinetics =  set_kinetics(p.network,'ms',parameters);

[t,c] = network_integrate(p.network,p.c_init,max(p.data.t));

internal   = find(p.network.external==0);

cc         = interp1(t,c(internal,:)',p.data.t)';
likelihood = exp( - 0.5 * sum(sum((cc-p.data.c).^2)) );

