function [res, es_result_list] = es_simulate_perturbation_sample(N,W,ind_ext,constraints,options,u_ratio,s_ratio,expansion,expansion_order,nrun)

%ES_SIMULATE_PERTURBATION_SAMPLE Predict the effects of a given differential profile 
%
% [res,es_result_list] = es_simulate_perturbation_sample(N,W,ind_ext,constraints,options,u_ratio,s_ratio,expansion,expansion_order,nrun)
%
% Compute the predicted effects of a given differential profile 
% for enzyme concentrations and external concentrations
% on concentrations and fluxes: do this for a number of model variants
% with sampled elasticities and present a statistics of the results
%
% N, W, ind_ext: structural information about the network
% x_ratio:       expression ratio column vector (or matrix, for several samples)
% nrun:          number of samples
% constraints, options: options for elasticity sampling
% s_ratio:       external concentration ratio vector
% u_ratio:       expression ratio vector
% expansion:     {'logarithmic','non-logarithmic'}
%                -> depending on the expansion type, the output refers to logarithmic 
%                   or non-logarithmic values; 
%
% Outputs:
% res:            contains statistics over predicted effects
% es_result_list: contains the results of all sampling runs

eval(default('expansion','''non-logarithmic''','expansion_order','2','nrun','1'));


% -------------------------------------------------
% iterate:
%   - sample elasticities 
%   - compute the metabolite and flux changes due to perturbation

for it = 1:nrun,
  display(sprintf('Monte Carlo run %d/%d',it,nrun));
  options.seed        = it;
  es_result_list{it}  = es_sample_model(N,W,ind_ext,constraints,options);
  r                   = es_simulate_perturbation(N,W,ind_ext,u_ratio,s_ratio,es_result_list{it},expansion,expansion_order);
  res.c(:,it)         = es_result_list{it}.c;
  res.v(:,it)         = es_result_list{it}.v;
  res.delta_c(:,:,it) = r.delta_c;
  res.delta_v(:,:,it) = r.delta_v;
  res.c_ratio(:,:,it) = r.c_ratio;
  res.v_ratio(:,:,it) = r.v_ratio;
end


% -------------------------------------------------
% statistics over the sampling results

flag_multiple_samples = double(size(u_ratio,2)>1);

res.c_mean = mean(res.c,2+flag_multiple_samples);
res.v_mean = mean(res.v,2+flag_multiple_samples);

if nrun > 1,
  res.c_std  = std(res.c,0,2+flag_multiple_samples);
  res.v_std  = std(res.v,0,2+flag_multiple_samples);
else 
  res.c_std  = nan * res.c_mean;
  res.v_std  = nan * res.v_mean;
end

res.delta_c_mean = squeeze(mean(res.delta_c,  3+flag_multiple_samples));
res.delta_c_std  = squeeze( std(res.delta_c,0,3+flag_multiple_samples));
res.delta_v_mean = squeeze(mean(res.delta_v,  3+flag_multiple_samples));
res.delta_v_std  = squeeze( std(res.delta_v,0,3+flag_multiple_samples));
res.c_ratio_mean = squeeze(mean(res.c_ratio,  3+flag_multiple_samples));
res.c_ratio_std  = squeeze( std(res.c_ratio,0,3+flag_multiple_samples));
res.v_ratio_mean = squeeze(mean(res.v_ratio,  3+flag_multiple_samples));
res.v_ratio_std  = squeeze( std(res.v_ratio,0,3+flag_multiple_samples));

