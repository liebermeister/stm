function [result, n_discarded] = stm_sample_stable_model(N, W, ind_ext, stm_constraints, stm_options)

% STM_SAMPLE_STABLE_MODEL - Run stm_sample_model until the steady state is stable
%
% [result, n_discarded] = stm_sample_stable_model(N, W, ind_ext, stm_constraints, stm_options)
%
% For the inputs, see stm_sample_model.
%
% stm_options.n_max_unstable: maximal number of unsuccessful attempts (default 100).
%                             If all of them yield unstable states, an error is thrown
%
% n_discarded  number of models discarded because of unstable states
%
% With a numerical seed in stm_options, the result is reproducible
% With stm_options.flag_control = 0, stability is not checked

stm_options = join_struct(struct('n_max_unstable',100),stm_options);

for n_discarded = 0:stm_options.n_max_unstable
  result = stm_sample_model(N, W, ind_ext, stm_constraints, stm_options);
  if ~isequal(result.control.stable,0), return; end
  stm_options.seed = nan; % otherwise, the same model would be sampled again
end

error('No stable model found: all %d attempts gave unstable states (limit: stm_options.n_max_unstable)',n_discarded+1);
