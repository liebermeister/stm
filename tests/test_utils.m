% Tests for helper functions (epistasis measures, statistics, options)

setup_test_path;
f_ref    = 1;
f_single = [0.8; 0.5];
f_double = [0.64 0.3; 0.3 0.25];

% Script-based tests: each section below is a separate test, and only the
% variables defined above are available in all sections

%% log_epistasis and simple_epistasis
e_log    = log_epistasis(f_ref, f_single, f_double);
e_simple = simple_epistasis(f_ref, f_single, f_double);
assert(max(max(abs(e_log - log(1 + e_simple)))) < 1e-12);
assert(abs(e_simple(1,2) - (0.3/0.4 - 1)) < 1e-12);

%% segre_epistasis: aggravating interaction, zero diagonal
e = segre_epistasis(f_ref, f_single, f_double);
assert(all(diag(e) == 0) && e(1,2) < 0);

%% multiple_testing_fdr
p = [0.001 0.002 0.5 0.9];
assert(isequal(multiple_testing_fdr(p, 0.05), logical([1 1 0 0])));

%% pvalue_from_countnumber
assert(pvalue_from_countnumber(10, 1:5) == 1/7);

%% randomize_synergies: symmetric, zero diagonal, same off-diagonal values
E = [0 1 2; 1 0 3; 2 3 0];
E_rnd = randomize_synergies(E);
assert(isequal(E_rnd, E_rnd') && all(diag(E_rnd) == 0));
assert(isequal(sort(E_rnd(triu(true(3),1))), [1 2 3]'));

%% es_constraints_set
network = network_construct([-1 0; 1 -1; 0 1], [1 1]', [1 3]');
[~, c] = es_default_options(network);
c = es_constraints_set(c, network, 'v_fix', network.actions(2), 3, 'v_sign', network.actions(1), 1, 'ext_signs', network.metabolites(3), -1);
assert(c.v_fix(2) == 3 && isnan(c.v_fix(1)) && c.v_sign(1) == 1 && c.ext_signs(3) == -1);

%% es_mcmc_prior_likelihood_sampling: samples follow the likelihood
rng(1);
prior      = @(p) struct('x', rand);
likelihood = @(s, p) exp(-(s.x - 0.8)^2 / 0.005);
[sample, ~, rate] = es_mcmc_prior_likelihood_sampling(prior, [], likelihood, [], 300, 20, {'x'});
x = cellfun(@(s) s.x, sample);
assert(abs(mean(x) - 0.8) < 0.05 && rate > 0 && rate < 1);
