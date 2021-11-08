function res = es_compare_ensembles(output1, output2, n_per, fdr, verbose)

% ES_COMPARE_ENSEMBLES - Compare an output function between two model ensembles
%
% res = es_compare_ensembles(output1, output2, n_per, fdr, verbose)
% 
% Significance test for model outputs obtained from multiple elasticity sampling
%
% Significant differences are computed from the outputs of a previous multiple 
% elasticity sampling for two model variants (see 'es_sample_ensemble').
% Both model variants must have the same number of metabolites and reactions.
%
% Inputs
%
%   output1  - Matrix or tensor from 1st model variant (from 'es_sample_ensemble')
%   output2  - Matrix or tensor from 2nd model variant (from 'es_sample_ensemble')
%   n_per    - Number of samples in permutation test
%   fdr      - False discovery rate
%
% Output: Data structure 'res' with fields:
%
%   res.mean_total               - Mean value for all samples (total mean)
%   res.mean_a                   - Mean value for 1st variant
%   res.mean_b                   - Mean value for 2nd variant
%   res.mean_delta               - Difference between values for both variants
%
%   res.p_value_mean_total       - p value of [mean value for all samples ~= total mean]
%   res.p_value_mean_a           - p value of [mean value for 1st variant ~= total mean]
%   res.p_value_mean_b           - p value of [mean value for 2nd variant ~= total mean]
%   res.p_value_mean_delta       - p value of [difference between values for both variants
%                                - ~= difference between values from shuffled variants]
%                                
%   res.mean_total_significant   - Significant elements for mean value for all samples
%   res.mean_a_significant       - Significant elements for mean value for 1st variant
%   res.mean_b_significant       - Significant elements for mean value for 2nd variant
%   res.mean_delta_significant   - Significant elements for difference between variants
%                                - Values: 1 (significantly high),  -1 (significantly low)
%
%   res.n_mean_total_significant - # significant elements for mean value for all samples
%   res.n_mean_a_significant     - # significant elements for mean value for 1st variant
%   res.n_mean_b_significant     - # significant elements for mean value for 2nd variant
%   res.n_mean_delta_significant - # significant elements for difference between variants


eval(default('n_per','100','fdr','0.05', 'verbose', '0'));

if norm(size(output1)-size(output1)), 
  error('function arguments output1 and output2 must have the same size'); 
end

switch length(size(output1)),
  
  case 2,
    D(:,:,1) = output1;
    D(:,:,2) = output2;
    res = influence_anova(D, n_per, fdr, [], verbose);
  
  case 3,
    D(:,:,:,1) = output1;
    D(:,:,:,1) = output2;
    res = interaction_anova(D, n_per, fdr, [], verbose);

end
