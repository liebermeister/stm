function p = interaction_anova(D,n_per,fdr,ind_intervention, verbose)

% function p = influence_anova(D, n_per, fdr, ind_intervention)
%
% Anova of a variable characterising pairs reactions 
%   (e.g., all second-order control coefficients for a single objective)
%
% Arguments:
% D:               data tensor (size nr x nr x n_sim x 2) where 
%                       nr: # reactions (corresponds to 'reaction' vector)
%                    n_sim: # number of MC samples
%                        2: two qualitative conditions to be compared (e.g. different fluxes)
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


eval(default('fdr','0.01','verbose','0'));

if exist('ind_intervention','var'), D = D(ind_intervention,ind_intervention,:,:); end

if verbose, display('Testing all pairs of interventions'); end

nr = size(D,1);
n_sim = size(D,3);
% to be sure, erase all values except for upper diagonal
for it = 1:size(D,3),
  Daa =  squeeze(D(:,:,it,1));
  Daa(find(tril(ones(nr)))) = nan;
  D(:,:,it,1) = Daa;
  Dbb =  squeeze(D(:,:,it,2));
  Dbb(find(tril(ones(nr)))) = nan;
  D(:,:,it,2) = Dbb;
end

Da = squeeze(D(:,:,:,1)); 
Db = squeeze(D(:,:,:,2)); 

finite_indices   = find(isfinite(D)); 
n_finite         = length(finite_indices);
finite_indices_a = find(isfinite(Da)); 
n_finite_a       = length(finite_indices_a);
finite_indices_b = find(isfinite(Db)); 
n_finite_b       = length(finite_indices_b);

p.mean_total   = nan * ones(nr);
p.p_value_mean_total       = nan * ones(nr);
p.mean_a = nan * ones(nr);
p.p_value_mean_a     = nan * ones(nr);
p.mean_b = nan * ones(nr);
p.p_value_mean_b     = nan * ones(nr);
p.mean_delta              = nan * ones(nr);
p.p_value_mean_delta      = nan * ones(nr);

for r1 = 1:nr-1,
  for r2 = r1+1:nr,
    
    if verbose, display(sprintf('Pair %d - %d',r1,r2)); end
    data_table = [squeeze(Da(r1,r2,:))'; squeeze(Db(r1,r2,:))'];
    
% mean interaction value
    
    p.mean_total(r1,r2) = mean(mean(data_table)');
    
% significance of the mean value (comparison with values drawn randomly
% from the entire data set)
    
    clear this_mean 
    
    for it_perm = 1:n_per, 
      rand_indices = randperm(n_finite);
      indices = finite_indices(rand_indices(1:2*n_sim));
      this_mean(it_perm) = mean(D(indices));
    end
    
    p.p_value_mean_total(r1,r2) =  pvalue_from_countnumber(p.mean_total(r1,r2),this_mean);

    % mean interaction value for condition A
    
    p.mean_a(r1,r2) = mean(data_table(1,:));
    
% significance of the mean value (comparison with values drawn randomly
% from the entire data set)
    
    clear this_mean_a 
    
    for it_perm = 1:n_per, 
      rand_indices = randperm(n_finite_a);
      indices = finite_indices(rand_indices(1:2*n_sim));
      this_mean(it_perm) = mean(Da(indices));
    end
    
    p.p_value_mean_a(r1,r2) =  pvalue_from_countnumber(p.mean_a(r1,r2),this_mean);

    
    % mean interaction value for condition B
    
    p.mean_b(r1,r2) = mean(data_table(2,:));
    
% significance of the mean value (comparison with values drawn randomly
% from the entire data set)
    
    clear this_mean_b 
    
    for it_perm = 1:n_per, 
      rand_indices = randperm(n_finite_b);
      indices = finite_indices(rand_indices(1:2*n_sim));
      this_mean(it_perm) = mean(Db(indices));
    end
    
    p.p_value_mean_b(r1,r2) = pvalue_from_countnumber(p.mean_b(r1,r2),this_mean);

% difference between two conditions
    p.mean_delta(r1,r2) = [1, -1] * mean(data_table,2);
    
% p-value for this difference to be larger than expected by chance
    
    clear this_dmean
    for it_perm = 1:n_per, 
      dummi = reshape(data_table(randperm(2*n_sim)),2,n_sim);
      this_dmean(it_perm) = [1, -1] * mean(dummi,2);
    end
    
    p.p_value_mean_delta(r1,r2) = pvalue_from_countnumber(p.mean_delta(r1,r2),this_dmean);
    
  end
end


% ---------------------------------------------------------------
% significance

p.mean_total_significant   = multiple_testing_fdr(p.p_value_mean_total,fdr)   - multiple_testing_fdr(1-p.p_value_mean_total  ,fdr);
p.mean_a_significant = multiple_testing_fdr(p.p_value_mean_a,fdr) - multiple_testing_fdr(1-p.p_value_mean_a,fdr);
p.mean_b_significant = multiple_testing_fdr(p.p_value_mean_b,fdr) - multiple_testing_fdr(1-p.p_value_mean_b,fdr);
p.mean_delta_significant  = multiple_testing_fdr(p.p_value_mean_delta,fdr)  - multiple_testing_fdr(1-p.p_value_mean_delta,fdr);

p.mean_total_significant(find(tril(ones(nr))))   = nan;
p.mean_a_significant(find(tril(ones(nr)))) = nan;
p.mean_b_significant(find(tril(ones(nr)))) = nan;
p.mean_delta_significant(find(tril(ones(nr))))  = nan;

p.n_mean_total_significant   = nansum(nansum(abs(p.mean_total_significant  )));
p.n_mean_a_significant       = nansum(nansum(abs(p.mean_a_significant)));
p.n_mean_b_significant       = nansum(nansum(abs(p.mean_b_significant)));
p.n_mean_delta_significant   = nansum(nansum(abs(p.mean_delta_significant )));
