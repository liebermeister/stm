function [alpha_A,alpha_I,alpha_M] = sample_saturation_parameters(N,W,ind_ext,es_options,es_constraints)

%SAMPLE_SATURATION_PARAMETERS - Sampling of saturation parameters
%
% [alpha_A,alpha_I,alpha_M] = sample_saturation_parameters(N,W,ind_ext,es_options)

if ~isfield(es_options,'seed'), es_options.seed = nan;     end
if ~isnan(es_options.seed), rand('state',es_options.seed); end

[Mplus, Mminus, Wplus, Wminus, nm, nr, N_int, ind_M, ind_Wp, ind_Wm, ind_Mp, ind_Mm] = make_structure_matrices(N,W,ind_ext,es_options.h);

if ~isfield(es_options,'set_alpha_to_half'), es_options.set_alpha_to_half = 0; end

empty = sparse(nr,nm); alpha_A = empty; alpha_I = empty; alpha_M = empty;

if es_options.set_alpha_to_half,
  alpha_M(ind_Mp) = 0.5;
  alpha_M(ind_Mm) = 0.5;  
  alpha_A(ind_Wp) = 0.5;
  alpha_I(ind_Wm) = 0.5;
  display(' Setting all alpha values = 0.5');
else,
  alpha_M(ind_Mp)  = rand(size(ind_Mp));
  alpha_M(ind_Mm)  = rand(size(ind_Mm));
  alpha_A(ind_Wp)  = rand(size(ind_Wp));
  alpha_I(ind_Wm)  = rand(size(ind_Wm));
  display(' Drawing alpha values uniformly from [0,1]');

  if es_options.set_alpha_nonuniform,
    %% THIS STILL HAS TO BE TESTED!!!!
    if size(es_constraints.alpha_M_mean),
      [Beta_M_alpha, Beta_M_beta] = Beta_mean_std_to_parameters(es_constraints.alpha_M_mean, es_constraints.alpha_M_std);
      rnd_M = betarnd(full(Beta_M_alpha), full(Beta_M_beta));
      alpha_M(ind_Mp(isfinite(rnd_M(ind_Mp)))) = rnd_M(ind_Mp(find(isfinite(rnd_M(ind_Mp)))));
      alpha_M(ind_Mm(isfinite(rnd_M(ind_Mm)))) = rnd_M(ind_Mm(find(isfinite(rnd_M(ind_Mm)))));
    end
    if size(es_constraints.alpha_A_mean),
      [Beta_A_alpha, Beta_A_beta] = Beta_mean_std_to_parameters(es_constraints.alpha_A_mean, es_constraints.alpha_A_std);
      rnd_A = betarnd(full(Beta_A_alpha), full(Beta_A_beta));
      alpha_A(ind_Wp(isfinite(rnd_A(ind_Wp)))) = rnd_A(ind_Wp(find(isfinite(rnd_A(ind_Wp)))));
    end
    if size(es_constraints.alpha_I_mean),
      [Beta_I_alpha, Beta_I_beta] = Beta_mean_std_to_parameters(es_constraints.alpha_I_mean, es_constraints.alpha_I_std);
      rnd_I = betarnd(full(Beta_I_alpha), full(Beta_I_beta));
      alpha_I(ind_Wm(isfinite(rnd_I(ind_Wm)))) = rnd_I(ind_Wm(find(isfinite(rnd_I(ind_Wm)))));
    end
    display(' Drawing some alpha values from Beta distribution, overriding uniform distribution');
  end
end

if sum(sum(isfinite([es_constraints.beta_M_fix(:),...
                 es_constraints.beta_A_fix(:),...
                 es_constraints.beta_I_fix(:)]))),
  display('Inserting given alpha values');
  alpha_M(isfinite(es_constraints.beta_M_fix)) = 1 - es_constraints.beta_M_fix(isfinite(es_constraints.beta_M_fix));
  alpha_A(isfinite(es_constraints.beta_A_fix)) = 1 - es_constraints.beta_A_fix(isfinite(es_constraints.beta_A_fix));
  alpha_I(isfinite(es_constraints.beta_I_fix)) = 1 - es_constraints.beta_I_fix(isfinite(es_constraints.beta_I_fix));
end
