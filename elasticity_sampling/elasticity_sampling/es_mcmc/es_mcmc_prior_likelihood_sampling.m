function [sample, likeli, acceptance_rate] =  es_mcmc_prior_likelihood_sampling(prior_function, prior_par, likelihood_function, likelihood_par,n_samples,n_relax,keep_fields)

%MCMC_PRIOR_LIKELIHOOD_SAMPLING - MCMC posterior sampling
%
% [sample, likeli, acceptance_rate] =  mcmc_prior_likelihood_sampling(prior_function, prior_par, likelihood_function, likelihood_par,n_samples,n_relax,keep_fields)
%
% mcmc posterior sampling using sampling from the prior 
% and acceptance/rejection according to the likelihood
%
% prior_function:      m-file that samples (once) from the prior
% prior_par            parameters for the prior sampling
% likelihood_function: m-file that computed the likelihood
% prior_par            parameters for the likelihood calculation
% n_samples:           # of samples
% n_relax:             # of samples before the actual sampling starts

sample = repmat({struct},n_relax+n_samples,1);

display(sprintf('\nInitial sample'));
old_sample = feval(prior_function, prior_par);
old_likeli = feval(likelihood_function, old_sample, likelihood_par);

accept_list = [];

for it = 1:n_relax+n_samples,
  display(sprintf('\nSample %d/%d',it-n_relax,n_samples));
  new_sample = feval(prior_function, prior_par);
  new_likeli = feval(likelihood_function, old_sample, likelihood_par);
  accept = 0;
  if new_likeli > old_likeli, 
    accept = 1;
  elseif rand < (new_likeli / old_likeli),
    accept = 1;
  end

  if accept, 
    old_sample = new_sample; 
    old_likeli = new_likeli; 
  end  

  accept_list(it) = accept;

  for it2 = 1:length(keep_fields),
    sample{it} = setfield(sample{it},keep_fields{it2},getfield(old_sample,keep_fields{it2}));
  end
  
  likeli(it) = old_likeli;

end

accept_list = accept_list(n_relax+1:end); 
sample = sample(n_relax+1:end);
likeli = likeli(n_relax+1:end);
acceptance_rate = sum(accept_list)/length(accept_list);
