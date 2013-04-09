function [alpha_A,alpha_I,alpha_M] = sample_saturation_parameters(N,W,ind_ext,options)

% [alpha_A,alpha_I,alpha_M] = sample_saturation_parameters(N,W,ind_ext,options)

if ~isfield(options,'seed'), options.seed = nan;     end
if ~isnan(options.seed), rand('state',options.seed); end

[Mplus, Mminus, Wplus, Wminus, nm, nr, N_int, ind_M, ind_Wp, ind_Wm, ind_Mp, ind_Mm] = make_structure_matrices(N,W,ind_ext,options.h);

if ~isfield(options,'set_alpha_to_half'), options.set_alpha_to_half = 0; end

empty = sparse(nr,nm); alpha_A = empty; alpha_I = empty; alpha_M = empty;

if options.set_alpha_to_half,
  alpha_A(ind_Wp) = 0.5;
  alpha_I(ind_Wm) = 0.5;
  alpha_M(ind_Mp) = 0.5;
  alpha_M(ind_Mm) = 0.5;  
  display(' Setting all alpha values = 0.5');
else,
  alpha_A(ind_Wp)  = rand(size(ind_Wp));
  alpha_I(ind_Wm)  = rand(size(ind_Wm));
  alpha_M(ind_Mp)  = rand(size(ind_Mp));
  alpha_M(ind_Mm)  = rand(size(ind_Mm));
  display(' Drawing alpha values uniformly from [0,1]');

  if options.set_alpha_nonuniform,
    %% THIS STILL HAS TO BE TESTED!!!!
    [Beta_A_alpha, Beta_A_beta] = Beta_mean_std_to_parameters(options.alpha_A_mean, options.alpha_A_std);
    [Beta_I_alpha, Beta_I_beta] = Beta_mean_std_to_parameters(options.alpha_A_mean, options.alpha_A_std);
    [Beta_M_alpha, Beta_M_beta] = Beta_mean_std_to_parameters(options.alpha_A_mean, options.alpha_A_std);
    rnd_A = betarnd(Beta_A_alpha, Beta_A_beta);
    rnd_I = betarnd(Beta_I_alpha, Beta_I_beta);
    rnd_M = betarnd(Beta_M_alpha, Beta_M_beta);
    alpha_A(isfinite(rnd_A(ind_Wp))) = rnd_A(find(isfinite(rnd_A(ind_Wp))));
    alpha_I(isfinite(rnd_I(ind_Wm))) = rnd_I(find(isfinite(rnd_A(ind_Wm))));
    alpha_M(isfinite(rnd_M(ind_Mp))) = rnd_M(find(isfinite(rnd_A(ind_Mp))));
    alpha_M(isfinite(rnd_M(ind_Mm))) = rnd_M(find(isfinite(rnd_A(ind_Mm))));
  end
end
