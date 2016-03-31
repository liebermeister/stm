function p = influence_anova(D, n_per, fdr, ind_intervention, verbose)

% function p = influence_anova(D, n_per, fdr, ind_intervention, verbose)
%
% Anova of a variable characterising single reactions 
%   (e.g., all first-order control coefficients for a single objective)
%
% Arguments:
% D:               data tensor (size nr x n_sim x 2) where 
%                       nr: # reactions (corresponds to 'reaction' vector)
%                    n_sim: # number of MC samples
%                        2: two qualitative conditions to be compared (e.g. different fluxes)
% n_per:            sample number in permutation tests
% fdr:              false discovery rate
% ind_intervention: reaction indices to be selected for the analysis (optional)
%
% The results are returned in a structure 'p' with fields
%
%   p.mean_total               Mean value for all samples (total mean)
%   p.mean_a                   Mean value for 1st variant
%   p.mean_b                   Mean value for 2nd variant
%   p.mean_delta               Difference between values for both variants
%
%   p.p_value_mean_total       p value of [mean value for all samples ~= total mean]
%   p.p_value_mean_a           p value of [mean value for 1st variant ~= total mean]
%   p.p_value_mean_b           p value of [mean value for 2nd variant ~= total mean]
%   p.p_value_mean_delta       p value of [difference between values for both variants
%                                           ~= difference between values from shuffled variants]
%
%   p.mean_total_significant   Significant elements for mean value for all samples
%   p.mean_a_significant       Significant elements for mean value for 1st variant
%   p.mean_b_significant       Significant elements for mean value for 2nd variant
%   p.mean_delta_significant   Significant elements for difference between values for both variants
%    Entries: 1 (for significantly high),  -1 (for significantly low)
%   
%   p.n_mean_total_significant # significant elements for mean value for all samples
%   p.n_mean_a_significant     # significant elements for mean value for 1st variant
%   p.n_mean_b_significant     # significant elements for mean value for 2nd variant
%   p.n_mean_delta_significant # significant elements for difference between values for both variants

eval(default('fdr','0.01','ind_intervention','[]','verbose','0'));

if length(ind_intervention), D = squeeze(D(ind_intervention,:,:)); end

if verbose,
  display('Testing all single interventions');
end

nr    = size(D,1);
n_sim = size(D,2);

Da = squeeze(D(:,:,1)); 
Db = squeeze(D(:,:,2)); 

finite_indices   = find(isfinite(D)); 
finite_indices_a = find(isfinite(Da)); 
finite_indices_b = find(isfinite(Db)); 
n_finite         = length(finite_indices);
n_finite_a       = length(finite_indices_a);
n_finite_b       = length(finite_indices_b);

p.mean_total         = nan * ones(nr,1);
p.mean_a             = nan * ones(nr,1);
p.mean_b             = nan * ones(nr,1);
p.mean_delta         = nan * ones(nr,1);
p.p_value_mean_total = nan * ones(nr,1);
p.p_value_mean_a     = nan * ones(nr,1);
p.p_value_mean_b     = nan * ones(nr,1);
p.p_value_mean_delta = nan * ones(nr,1);


for r1 = 1:nr,
    
    if verbose, display(sprintf('Reaction %d',r1)); end
    data_table = [Da(r1,:); Db(r1,:)];
    
    %% total mean value
    %% p-value: background distribution sampled from condition A and B
    %% this test should typically lead to INSIGNIFICANT results!
    
    p.mean_total(r1) = mean(mean(data_table));
    this_mean = [];
    for it_perm = 1:n_per, 
      rand_indices = randperm(n_finite);
      indices      = finite_indices(rand_indices(1:2*n_sim));
      this_mean(it_perm) = mean(D(indices));
    end
    p.p_value_mean_total(r1) = pvalue_from_countnumber(p.mean_total(r1),this_mean);

    %% mean value for condition A
    %% p-value: background distribution sampled from condition A and B

    p.mean_a(r1) = mean(data_table(1,:));    
    this_mean_a = [];
    for it_perm = 1:n_per, 
      rand_indices = randperm(n_finite_a);
      indices = finite_indices(rand_indices(1:2*n_sim));
      this_mean(it_perm) = mean(Da(indices));
    end
    p.p_value_mean_a(r1) =  pvalue_from_countnumber(p.mean_a(r1),this_mean);
    
    %% mean value for condition B
    %% p-value: background distribution sampled from condition A and B
    
    p.mean_b(r1) = mean(data_table(2,:));    
    this_mean_b = [];
    for it_perm = 1:n_per, 
      rand_indices = randperm(n_finite_b);
      indices = finite_indices(rand_indices(1:2*n_sim));
      this_mean(it_perm) = mean(Db(indices));
    end
    p.p_value_mean_b(r1) =  pvalue_from_countnumber(p.mean_b(r1),this_mean);
    
    %% difference between two conditions
    %% p-value: background distribution sampled from condition A and B
    
    p.mean_delta(r1) = [1, -1] * mean(data_table,2);    
    this_mean_delta = [];
    for it_perm = 1:n_per, 
      dummi = reshape(data_table(randperm(2*n_sim)),2,n_sim);
      this_mean_delta(it_perm) = [1, -1] * mean(dummi,2);
    end
    p.p_value_mean_delta(r1) = pvalue_from_countnumber(p.mean_delta(r1),this_mean_delta);
    
  end

  
% ---------------------------------------------------------------
% determine significant elements

p.mean_total_significant = multiple_testing_fdr(p.p_value_mean_total,fdr) ...
    - multiple_testing_fdr(1-p.p_value_mean_total  ,fdr);

p.mean_a_significant = multiple_testing_fdr(p.p_value_mean_a,fdr) ...
    - multiple_testing_fdr(1-p.p_value_mean_a,fdr);,

p.mean_b_significant = multiple_testing_fdr(p.p_value_mean_b,fdr) ...
    - multiple_testing_fdr(1-p.p_value_mean_b,fdr);

p.mean_delta_significant = multiple_testing_fdr(p.p_value_mean_delta,fdr) ...
    - multiple_testing_fdr(1-p.p_value_mean_delta,fdr);


% ---------------------------------------------------------------
% count significant elements

p.n_mean_total_significant = sum(sum(abs(p.mean_total_significant)));
p.n_mean_a_significant     = sum(sum(abs(p.mean_a_significant)));
p.n_mean_b_significant     = sum(sum(abs(p.mean_b_significant)));
p.n_mean_delta_significant = sum(sum(abs(p.mean_delta_significant )));
