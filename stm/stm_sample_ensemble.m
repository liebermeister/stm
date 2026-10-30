function [output, output_list] = stm_sample_ensemble(N, W, ind_ext, stm_constraints, stm_options, output_function, select_function, score_function, function_args)

% STM_SAMPLE_ENSEMBLE - Create model ensemble and collect values for an output function
%
% [output, output_list] = stm_sample_ensemble(N, W, ind_ext, stm_constraints, stm_options, output_function, select_function, score_function, function_args)
%
% Run elasticity sampling repeatedly, compute every time an output function
% (any kind of matlab variable), and return a list of the output values.
%
% This function also allows for posterior sampling (Metropolis-Monte Carlo algorithm)
%
% There is no direct way to run the sampling in parallel for two scenarios
% with the same random values used for saturation constants in both ensembles
%
% Inputs (nm: # metabolites; nr: # reactions)
%   N         - Stoichiometric matrix (nm x nr)
%   W         - Allosteric regulation matrix (nr x nm)
%   ind_ext   - indices of external metabolites
%
%  For the inputs stm_constraints and stm_options, see stm_default_options
%
%  output_function - Output function for which response coefficients are sampled (*1)
%  select_function - Function for filtering the sampled models (*2)
%  score_function  - Function for posterior probability sampling (*3)
%  function_args   - Arguments for output_function, select_function, score_function
%
% Outputs
%  output_list     - Results as a list of matrices
%  output          - Same results as a tensor
%
% Comments:
%
% (1) The output function is called as output_function(result, function_args).
%     If its value is numeric, 'output' contains all values as a column vector
%     (for scalars), matrix (for row or column vectors), or tensor (for matrices).
%
% (2) If a 'select_function' is provided, only samples that yield a non-zero value
%     of this function will be considered.
%
% (3)  If a score function (a log likelihood) is provided, it is evaluated
%      every time for each sample and the values are used for Metropolis sampling
%      instead of simple repeated sampling.


eval(default('select_function','[]','score_function','[]','function_args','[]'));

if stm_options.set_alpha_to_half
  warning('Attempted ensemble sampling with fixed saturation values');
end

% -----------------------------------------------
% run sampling algorithm

old_score = -inf;

stm_options.seed = nan;

for it = 1:stm_options.n_samples
  display(sprintf('\nSample %d/%d',it,stm_options.n_samples));
  my_result = stm_sample_model(N,W,ind_ext,stm_constraints,stm_options);
  if ~isempty(select_function)
    accept = feval(select_function,my_result,function_args);
    while ~accept
      display(' Invalid sample; drawing another sample');
      my_result = stm_sample_model(N,W,ind_ext,stm_constraints,stm_options);
      accept = feval(select_function,my_result,function_args);
    end
  end
  my_output = feval(output_function,my_result,function_args);
  if isempty(score_function)
    output_list{it} = my_output;
  else
    my_score = feval(score_function,my_result,function_args);
    if [my_score > old_score] + [rand < exp(my_score - old_score)]
      output_list{it} = my_output;
      old_score = my_score;
    else
      output_list{it} = output_list{it-1};
    end
  end
end

% -----------------------------------------------
% put the output list into a vector / matrix / tensor

output = [];

if isscalar(output_list{1})

  for it = 1:length(output_list)
    output(it,:) = output_list{it};
  end

else

  if isvector(output_list{1})
    for it = 1:length(output_list)
      output(:,it) = output_list{it};
    end

  else % tensor

    if isnumeric(output_list{1})
      for it = 1:length(output_list)
        output(:,:,it) = output_list{it};
      end
    end
  end

end
