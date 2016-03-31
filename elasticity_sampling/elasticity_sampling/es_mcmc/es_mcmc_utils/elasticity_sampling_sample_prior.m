function result = elasticity_sampling_sample_prior(p)

% result = elasticity_sampling_sample_prior(p)

result = es_sample_model(p.N,p.W,p.ext_ind,p.constraints,p.options);

