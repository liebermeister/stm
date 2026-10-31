% Regression tests for errors found in the code review of October 2026 (see CHANGELOG.md)

setup_test_path;
[N, W, ind_ext, v, stm_options, stm_constraints] = example_model;
N_chain   = [-1 0; 1 -1; 0 1];
W_chain   = zeros(2,3);
ext_chain = [1 3]';
N_loop    = [-1 0 0 0 0; 1 -1 0 -1 0; 0 1 -1 0 0; 0 0 1 1 -1; 0 0 0 0 1]; % rank-deficient
N_branch  = [-1 0 0 0; 1 -1 -1 0; 0 1 0 -1; 0 0 1 0; 0 0 0 1];
N_fb      = [-1 0; 1 -1; 0 1];  % product activates its own synthesis
W_fb      = [0 2 0; 0 0 0];

%% rank-deficient stoichiometric matrix: thermodynamic forces are consistent
[o, c] = stm_default_options(N_loop); o.flag_second_order = 0; o.seed = 1;
c.v_fix = [2 1 1 1 2]';
r = stm_sample_model(N_loop, zeros(5,5), [1 5]', c, o);
assert(max(abs(N_loop' * r.mu + r.A)) < 1e-8);
assert(all(sign(r.A) == sign(c.v_fix)));

%% cooperativity h: rate law and elasticities agree with the stated thermodynamics
[o, c] = stm_default_options(N_chain); o.flag_second_order = 0; o.seed = 1; o.h = [2 0.5]';
for law = {'cs','ms'}
  o.kinetic_law = law{1};
  c.v_fix = [1 1]';
  r = stm_sample_model(N_chain, W_chain, ext_chain, c, o);
  k = r.kinetics;
  f = @(x) modular_velocities(k.type, N_chain, W_chain, ext_chain, k.u, x, k.KA, k.KI, k.KM, k.KV, k.Keq, k.h);
  assert(max(abs(f(r.c) - r.v)) < 1e-8);
  e = 1e-6; Ef = zeros(2,3);
  for i = 1:3
    d = zeros(3,1); d(i) = e * r.c(i);
    Ef(:,i) = (f(r.c + d) - f(r.c - d)) / (2 * e * r.c(i));
  end
  assert(max(abs(Ef(:) - full(r.elasticities.un_E_c(:)))) < 1e-6);
  c_eq = exp(pinv(N_chain') * log(k.Keq));
  assert(max(abs(f(c_eq))) < 1e-8);
end

%% large network with a conserved moiety: rank of N is computed correctly
nr = 110;
N_m = zeros(nr+3, nr);
for i = 1:nr, N_m(i,i) = -1; N_m(i+1,i) = 1; end
N_m(nr+2,5) = -1; N_m(nr+3,5) = 1; N_m(nr+2,6) = 1; N_m(nr+3,6) = -1;
ext = [1 nr+1]';
[o, c] = stm_default_options(N_m); o.flag_second_order = 0; o.seed = 1; o.verbose = 0;
c.v_fix = ones(nr,1); c.dmu_fix = -5 * ones(nr,1);
r = stm_sample_model(N_m, zeros(nr,nr+3), ext, c, o);
CS = full(r.control.CS); CJ = full(r.control.CJ);
assert(all(isfinite(CS(:))) && max(abs(CS(:))) < 1e6);
assert(max(abs(CJ * ones(nr,1) - 1)) < 1e-6);
assert(max(abs(CS(nr+2,:) + CS(nr+3,:))) < 1e-6); % moiety is conserved

%% predefined mu_eqconstraint: the sign of the thermodynamic forces is kept
[o, c] = stm_default_options(N_chain); o.flag_second_order = 0; o.seed = 1;
c.v_fix = [1 1]'; c.log_c_fix = zeros(3,1); c.dmu_fix = [-5 -7]';
c.mu_eqconstraint = struct('matrix',[1 0 0]','vector',0);
[~,~,~,~,~,mu,~,~,~,A] = stm_sample_steady_state(N_chain, W_chain, ext_chain, c, o);
assert(max(abs(N_chain' * mu + A)) < 1e-8);
assert(abs(mu(1)) < 1e-8);

%% very small mass action ratios are not clipped
[o, c] = stm_default_options(N_chain); o.flag_second_order = 0; o.seed = 1;
c.v_fix = [1 1]'; c.log_c_fix = [0 -15 -30]'; c.dmu_fix = [-5 -7]';
[~,~,~,~,Keq,mu,~,~,~,A] = stm_sample_steady_state(N_chain, W_chain, ext_chain, c, o);
assert(max(abs(A - RT * (log(Keq) - N_chain' * c.log_c_fix))) < 1e-8);
assert(max(abs(A - [5 7]')) < 1e-8);

%% a seed makes results reproducible, but model variables are not built from identical random numbers
[o, c] = stm_default_options(N); o.flag_second_order = 0; o.seed = 3;
c.v_fix = v; c.log_c_std = ones(6,1); c.log_u_std = ones(5,1);
r1 = stm_sample_model(N, W, ind_ext, c, o);
r2 = stm_sample_model(N, W, ind_ext, c, o);
assert(isequal(r1.KV, r2.KV));
assert(abs(corr(log(r1.c(1:5)), log(r1.u))) < 0.99);

%% solve_quadratic2 keeps the state of the random number generator
rng(5); s1 = rng;
solve_quadratic2(eye(2), 1, [1 1], zeros(1,2,2));
s2 = rng;
assert(isequal(s1.State, s2.State));

%% solve_quadratic2 finds the minimum
A = eye(3); B = [1 1 1; 0 0 1]; c = [1 2]';
D = zeros(2,3,3); D(1,:,:) = [1 .3 0; .3 2 .1; 0 .1 1]; D(2,:,:) = [.5 0 .2; 0 1 0; .2 0 3];
x = solve_quadratic2(A, c, B, D);
g = [x'*squeeze(D(1,:,:))*x; x'*squeeze(D(2,:,:))*x]/2 + B*x - c;
assert(norm(g) < 1e-6);
assert(abs(.5 * x' * A * x - 0.39673) < 1e-4);

%% stm_simulate_perturbation_sample: mean values have the size of one simulation
[~, ~, ~, ~, o, c] = example_model;
o.seed = nan; rng(1);
res = stm_simulate_perturbation_sample(N, W, ind_ext, c, o, [1.1 1; 1 1.1; 1 1; 1 1; 1 1], ones(4,2), 'non-logarithmic', 1, 5);
assert(isequal(size(res.delta_v_mean), [5 2]));
assert(size(res.c_mean,1) == 6);

%% check_hyperantagonistic: solutions satisfy the definition
rng(2);
for i = 1:300
  p = randn(1,5);
  [~, s, sol] = evalc('check_hyperantagonistic(p(1),p(2),p(3),p(4),p(5))');
  if s
    f = @(a,b) p(1)*a + p(2)*b + .5*(p(3)*a^2 + 2*p(4)*a*b + p(5)*b^2);
    a = sol(1); b = sol(2);
    assert(f(-a,0) < 0 && f(0,-b) < 0 && f(-a,-b) > f(-a,0) && f(-a,-b) > f(0,-b));
  end
end

%% pvalue_from_countnumber: missing value gives a missing p value
assert(isnan(pvalue_from_countnumber(nan, 1:5)));

%% fdr_two_sided: signs, and nan values are not significant
p = [1e-6 1-1e-6 0.5 nan];
assert(isequal(fdr_two_sided(p, 0.05), [1 -1 0 0]));

%% stm_compare_ensembles: outputs of scalar functions
rng(1);
r = stm_compare_ensembles(randn(40,1), randn(40,1) + 5, 200, 0.05);
assert(isscalar(r.p_value_mean_delta));

%% segre_epistasis: a lethal mutation gives no nan values
e = segre_epistasis(1, [1; 0; 0.8], [1 0 .8; 0 0 0; .8 0 .64]);
assert(~any(isnan(e(:))));

%% multi-reactions (more than 6 metabolites): normalised elasticities, also with second-order terms
N_multi = [-ones(6,1) zeros(6,1); 1 -1; 0 1];
[o, c] = stm_default_options(N_multi); o.seed = 1; o.flag_second_order = 1;
o.limit_elasticities_in_multireactions = 1;
c.v_fix = [1 1]';
r = stm_sample_model(N_multi, zeros(2,8), [1:6 8]', c, o);
assert(abs(sum(abs(full(r.elasticities.sc_E_c(1,:)))) - 1) < 1e-8);
assert(norm(double(full(r.elasticities.un_E_cc(1,:,:)))) == 0);
assert(all(isfinite(full(r.elasticities.un_E_c(:)))));

%% inactive reaction with a thermodynamic force, no_equilibrium = 1: enzyme is switched off
for no_equilibrium = [1 0]
  [o, c] = stm_default_options(N_branch); o.flag_second_order = 0; o.seed = 1;
  o.no_equilibrium = no_equilibrium;
  c.v_fix = [1 1 0 1]'; c.log_c_fix = zeros(5,1); c.dmu_fix = [-5 -5 -3 -5]';
  r = stm_sample_model(N_branch, zeros(4,5), [1 4 5]', c, o);
  if no_equilibrium
    assert(r.u(3) == 0 && r.KV(3) > 0);
  else
    assert(r.u(3) > 0 && r.KV(3) == 0);
  end
  assert(r.v(3) == 0);
  % scaled coefficients of the inactive reaction itself are undefined
  assert(~any(isnan(reshape(r.control.CJ_sc([1 2 4],:),1,[]))));
  assert(~any(isnan(reshape(r.control.RJu_sc([1 2 4],:),1,[]))));
  assert(all(isfinite(full(r.control.CJ(:)))));
end

%% reaction at equilibrium: elasticities agree with the rate law
[o, c] = stm_default_options(N_branch); o.flag_second_order = 0; o.seed = 1;
c.v_fix = [1 1 0 1]'; c.log_c_fix = zeros(5,1); c.dmu_fix = [-5 -5 0 -5]';
r = stm_sample_model(N_branch, zeros(4,5), [1 4 5]', c, o);
k = r.kinetics;
f = @(x) modular_velocities(k.type, N_branch, zeros(4,5), [1 4 5]', k.u, x, k.KA, k.KI, k.KM, k.KV, k.Keq, k.h);
e = 1e-6; d = [0 e*r.c(2) 0 0 0]';
fd = (f(r.c + d) - f(r.c - d)) / (2 * e * r.c(2));
assert(abs(fd(3)) > 1e-3);
assert(abs(fd(3) - full(r.elasticities.un_E_c(3,2))) < 1e-6);

%% spectral_density_to_smooth_variance: Ornstein-Uhlenbeck process, variance 1/(2a)
a = 1; w = logspace(-4,4,4000);
V = spectral_density_to_smooth_variance(w, 1./(w.^2 + a^2), 1e-6);
assert(abs(V - 1/(2*a)) < 1e-3);

%% energetics_control_and_noise uses the given concentrations
assert(contains(fileread(which('energetics_control_and_noise')), 'log_c_fix     = log(c)'));

%% thermodynamic forces: uniform sampling by default, vertices as an alternative
[o, c] = stm_default_options(N_chain); o.seed = nan; o.flag_second_order = 0;
c.v_fix = [1 1]'; c.log_c_fix = zeros(3,1); rng(1);
n = 300; A_all = zeros(2,n); A_ext = zeros(2,n);
o2 = o; o2.mu_sampling_method = 'extreme_points';
for i = 1:n
  [~,~,~,~,~,~,~,~,~,A_all(:,i)]  = stm_sample_steady_state(N_chain, W_chain, ext_chain, c, o);
  [~,~,~,~,~,~,~,~,~,A_ext(:,i)] = stm_sample_steady_state(N_chain, W_chain, ext_chain, c, o2);
end
assert(all(A_all(:) > 0) && all(A_all(:) <= 50));
assert(max(abs(mean(A_all,2) - 25)) < 3);
assert(min(std(A_all,0,2)) > 10);        % uniform distribution: 50/sqrt(12) = 14
assert(numel(unique(round(A_all(1,:),6))) > 250);
assert(numel(unique(round(A_ext(1,:),6))) < 10); % finite number of vertices

%% compute_modular_response_second agrees with the full second-order response coefficients
[o, c] = stm_default_options(N); o.flag_second_order = 1; o.seed = 2;
o.set_alpha_to_half = 0; o.kinetic_law = 'cs'; c.v_fix = v;
r = stm_sample_model(N, W, ind_ext, c, o);
s = r.saturation;
[R1, R2] = compute_modular_response_second([], [1 0 0 0 0]', o.kinetic_law, N, W, ind_ext, s.alpha_A, s.alpha_I, s.alpha_M, r.v, r.A, r.u, r.c, o.h, r.v_plus, r.v_minus, r.control.CS, r.control.CJ);
RJuu = double(full(r.control.RJuu_sc));
assert(max(abs(full(R1(:)) - full(reshape(r.control.RJu_sc(1,:),[],1)))) < 1e-10);
assert(max(max(abs(full(R2) - squeeze(RJuu(1,:,:))))) < 1e-10);

%% second-order response coefficients in a small network (dense tensors)
[o, c] = stm_default_options(N_chain); o.flag_second_order = 1; o.seed = 1; c.v_fix = [1 1]';
r = stm_sample_model(N_chain, W_chain, ext_chain, c, o);
assert(size(r.control.RJpp_un,1) == 2);

%% unstable states: nan control coefficients, stable models are drawn on request
[o, c] = stm_default_options(N_fb); o.set_alpha_to_half = 0; o.seed = nan; rng(1);
c.v_fix = [1 1]';
stable = zeros(100,1);
for flag_second_order = [0 1]
  o.flag_second_order = flag_second_order;
  for i = 1:100
    r = stm_sample_model(N_fb, W_fb, [1 3]', c, o);
    stable(i) = r.control.stable;
    if ~stable(i)
      assert(all(isnan(r.control.CS(:))) && all(isnan(r.control.CJ(:))) && all(isnan(r.control.RJp_sc(:))));
    end
  end
  assert(any(stable == 0) && any(stable == 1));
end
o.flag_second_order = 0;
rng(1);
n_discarded = 0;
for i = 1:30
  [r, n] = stm_sample_stable_model(N_fb, W_fb, [1 3]', c, o);
  n_discarded = n_discarded + n;
  assert(r.control.stable == 1);
end
assert(n_discarded > 0);

%% unstable states: a numerical seed does not cause an endless loop, and is reproducible
[o, c] = stm_default_options(N_fb); o.set_alpha_to_half = 0; o.flag_second_order = 0; c.v_fix = [1 1]'; o.seed = 1;
[r1, n1] = stm_sample_stable_model(N_fb, W_fb, [1 3]', c, o);
[r2, n2] = stm_sample_stable_model(N_fb, W_fb, [1 3]', c, o);
assert(n1 == n2 && isequal(r1.KV, r2.KV));

%% stm_sample_ensemble contains only stable models
[o, c] = stm_default_options(N_fb); o.set_alpha_to_half = 0; o.flag_second_order = 0; c.v_fix = [1 1]'; o.seed = nan; o.n_samples = 20; rng(1);
out = stm_sample_ensemble(N_fb, W_fb, [1 3]', c, o, @(r,a) double(r.control.stable));
assert(all(out == 1));

%% flag_control = 0: no control analysis
[o, c] = stm_default_options(N_fb); o.set_alpha_to_half = 0; c.v_fix = [1 1]'; o.seed = 1; o.flag_control = 0;
r = stm_sample_model(N_fb, W_fb, [1 3]', c, o);
assert(isempty(r.control.stable) && isfield(r,'elasticities'));

%% default limit for thermodynamic forces
[~, c] = stm_default_options(N_chain);
assert(c.dmu_limit == 50);

%% forces above dmu_limit are accepted if they follow from given data
[o, c] = stm_default_options(N_chain); o.flag_second_order = 0; o.seed = 1;
c.v_fix = [1 1]'; c.log_c_fix = zeros(3,1); c.Keq_fix = exp([60 5]' / RT);
[~,~,~,~,~,~,~,~,~,A] = stm_sample_steady_state(N_chain, W_chain, ext_chain, c, o);
assert(abs(A(1) - 60) < 1e-8);

%% stm_sample_stable_model: error after the maximal number of unsuccessful attempts
[o, c] = stm_default_options(N_fb); o.set_alpha_to_half = 0; o.flag_second_order = 0; c.v_fix = [1 1]'; o.seed = nan;
assert(o.n_max_unstable == 100);

o.n_max_unstable = 0;
rng(1); n_error_0 = 0;
for i = 1:100
  try
    [~, n] = stm_sample_stable_model(N_fb, W_fb, [1 3]', c, o);
    assert(n == 0);
  catch e
    assert(contains(e.message, 'n_max_unstable'));
    n_error_0 = n_error_0 + 1;
  end
end
assert(n_error_0 > 0); % with no redraws allowed, every unstable sample gives an error

%% fdr_check_permutations: warning if too few permutations
lastwarn('');
assert(fdr_check_permutations(rand(100,1), 100, 0.05) == 40);
assert(contains(lastwarn, 'too few'));
lastwarn('');
assert(fdr_check_permutations([rand(100,1); nan], 3998, 0.05) == 1);
assert(isempty(lastwarn));
