% Tests for elasticity sampling (steady state, elasticities, control coefficients)

setup_test_path;
[N, W, ind_ext, v, stm_options, stm_constraints] = example_model;
[nm, nr] = size(N);
internal = setdiff(1:nm, ind_ext);
K = null(N(internal,:));
result = stm_sample_model(N, W, ind_ext, stm_constraints, stm_options);
CJ = full(result.control.CJ);
CS = full(result.control.CS);

% Script-based tests: each section below is a separate test, and only the
% variables defined above are available in all sections

%% stm_default_options: fields used by the sampling functions
[o, c] = stm_default_options(N);
assert(isfield(c, 'ext_signs') && isequal(size(c.ext_signs), [nm 1]));
assert(isfield(o, 'limit_elasticities_in_multireactions'));

%% stm_sample_model: model reproduces the given fluxes, thermodynamics, and MCA theorems
assert(max(abs(result.v - v)) < 1e-8);
assert(all(sign(result.A) == sign(v)));
assert(max(abs(log(result.Keq) - N' * log(result.c0))) < 1e-8);
assert(max(abs(result.A - RT * log(result.Keq ./ result.Kma))) < 1e-8);
assert(max(max(abs(CJ * K - K))) < 1e-8);
assert(max(max(abs(CS * K))) < 1e-8);

%% stm_sample_model: rate laws in result.kinetics reproduce the fluxes
k = result.kinetics;
v_kin = modular_velocities(k.type, N, W, ind_ext, k.u, k.c, k.KA, k.KI, k.KM, k.KV, k.Keq, k.h);
assert(max(abs(v_kin - v)) < 1e-8);

%% stm_sample_steady_state: given equilibrium constants and concentrations
es_c = stm_constraints;
es_c.log_c_fix = log([1 1 0.5 0.2 0.1 0.1]');
es_c.Keq_fix   = exp(N' * es_c.log_c_fix + 2 * sign(v));
[c0, c, u, J, Keq, mu, Kc, Kma, zeta, A] = stm_sample_steady_state(N, W, ind_ext, es_c, stm_options);
assert(max(abs(A - 2 * RT * sign(v))) < 1e-8);
assert(max(abs(N' * mu + A)) < 1e-8);

%% stm_sample_model: second-order response coefficients
es_o = stm_options;
es_o.flag_second_order = 1;
result2 = stm_sample_model(N, W, ind_ext, stm_constraints, es_o);
assert(isfield(result2.control, 'RSpp_sc') && isfield(result2.control, 'RJp_un'));
assert(max(max(abs(full(result2.control.CJ) - CJ))) < 1e-8);

%% stm_sample_ensemble: shape of the output
es_o = stm_options;
es_o.n_samples = 3;
es_o.set_alpha_to_half = 0;
es_o.verbose = 0;
out = stm_sample_ensemble(N, W, ind_ext, stm_constraints, es_o, @(r, a) sum(r.control.CJ(:).^2));
assert(isequal(size(out), [3 1]));
out = stm_sample_ensemble(N, W, ind_ext, stm_constraints, es_o, @(r, a) r.KV);
assert(isequal(size(out), [nr 3]));
out = stm_sample_ensemble(N, W, ind_ext, stm_constraints, es_o, @(r, a) full(r.control.CJ));
assert(isequal(size(out), [nr nr 3]));

%% stm_compare_ensembles: differences between ensembles (vector outputs)
rng(1);
out1 = randn(4, 30);
out2 = out1 + 0.01 * randn(4, 30);
out2(1,:) = out2(1,:) + 10;
res = stm_compare_ensembles(out1, out2, 200, 0.05);
assert(abs(res.mean_delta(1) + 10) < 2);
assert(res.mean_delta_significant(1) == -1);
assert(all(res.mean_delta_significant(2:4) == 0));

%% stm_compare_ensembles: matrix outputs
out1 = randn(3, 3, 20);
out2 = randn(3, 3, 20);
out2(1,2,:) = out2(1,2,:) + 10;
res = stm_compare_ensembles(out1, out2, 100, 0.05);
assert(isequal(size(res.mean_delta), [3 3]));
assert(abs(res.mean_delta(1,2) + 10) < 2);
