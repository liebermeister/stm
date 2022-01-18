function [ris, es_result_list] = es_discover_perturbation_sample(N,W,ext_ind,constraints,options,c_ratio,v_ratio,expansion,expansion_order,nrun,u_sigma_prior,s_known_sigma_prior,s_unknown_sigma_prior)

%ES_DISCOVER_PERTURBATION_SAMPLE
%
% [ris, es_result_list] = es_discover_perturbation_sample(N,W,ext_ind,constraints,...
%  options,c_ratio,v_ratio,expansion,expansion_order,nrun,u_sigma_prior,s_sigma_prior)
%
% Compute the most likely perturbation of enzyme expression (u) and external metabolites (s)
% that would give rise to a given differential profile of balanced concentrations c and fluxes j
% do this for a number of model variants
% with sampled elasticities and present a statistics of the results
%
% N, W, ext_ind: structural information about the network
% x_ratio:       expression ratio column vector (or matrix, for several samples)
% nrun:          number of samples
% constraints, options: options for elasticity sampling
% expansion:     {'logarithmic','non-logarithmic'}
%                -> depending on the expansion type, the output refers to logarithmic 
%                   or non-logarithmic values; 
% expansion_order:       1 or 2 for first- or second-expansion_order expansion
%
% Outputs:
% ris:            contains statistics over predicted effects
% es_result_list: contains the results of all sampling runs

eval(default('expansion','''non-logarithmic''','expansion_order','2','nrun','1'));

[nm,nr] = size(N);

options.set_alpha_to_half = 0;

% -------------------------------------------------
% iterate:
%   - sample elasticities 
%   - compute the metabolite and flux changes due to perturbation

ris.u_ratio = [];
ris.s_ratio = [];
ris.c_ratio = [];
ris.v_ratio = [];

for it = 1:nrun,
  
  display(sprintf('Monte Carlo run %d',it));

  this_es_result      = es_sample_model(N,W,ext_ind,constraints,options);

  [u_ratio, s_ratio, pred_effect] = es_discover_perturbation(...
      N,W,ext_ind,c_ratio,v_ratio,this_es_result,expansion,expansion_order,...
      u_sigma_prior,s_known_sigma_prior,s_unknown_sigma_prior);

  ris.u_ratio(:,it)   = u_ratio;
  ris.s_ratio(:,it)   = s_ratio;
  ris.c_ratio(:,it)   = pred_effect.c_ratio;
  ris.v_ratio(:,it)   = pred_effect.v_ratio;
  if nargout>1,   es_result_list{it} = this_es_result; end

end


% -------------------------------------------------
% statistics over the sampling results

ris.u_ratio_mean = mean(ris.u_ratio,2);
ris.s_ratio_mean = mean(ris.s_ratio,2);
ris.c_ratio_mean = mean(ris.c_ratio,2);
ris.v_ratio_mean = mean(ris.v_ratio,2);

if nrun > 1,
  ris.u_ratio_std  = std(ris.u_ratio,0,2);
  ris.s_ratio_std  = std(ris.s_ratio,0,2);
  ris.c_ratio_std  = std(ris.c_ratio,0,2);
  ris.v_ratio_std  = std(ris.v_ratio,0,2);
else 
  ris.u_ratio_std  = nan*ris.u_ratio_mean;
  ris.s_ratio_std  = nan*ris.s_ratio_mean;
  ris.c_ratio_std  = nan*ris.c_ratio_mean;
  ris.v_ratio_std  = nan*ris.v_ratio_mean;
end


% -------------------------------------------------
% probabilities for being the most strongly downregulated enzyme

[dum,ranks] = sort(ris.u_ratio,1);

if nrun>0,
for it = 1:nr, ris.prob_to_be_best(it) = sum(ranks(1,:)==it)/nrun; end
end