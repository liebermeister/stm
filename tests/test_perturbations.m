% Tests for the prediction of perturbation effects

setup_test_path;
[N, W, ind_ext, v, es_options, es_constraints] = example_model;
[nm, nr] = size(N);
es_options.flag_second_order = 1;
result = es_sample_model(N, W, ind_ext, es_constraints, es_options);
s_ratio = ones(length(ind_ext), 1);
u_ratio = [1.01 1 1 0.99 1]';
r1 = es_simulate_perturbation(N, W, ind_ext, u_ratio, s_ratio, result, 'logarithmic', 1);
r2 = es_simulate_perturbation(N, W, ind_ext, u_ratio, s_ratio, result, 'logarithmic', 2);
r3 = es_simulate_perturbation(N, W, ind_ext, u_ratio, s_ratio, result, 'non-logarithmic', 1);
r4 = es_simulate_perturbation(N, W, ind_ext, u_ratio, s_ratio, result, 'non-logarithmic', 2);
A = eye(3); B = [1 1 1; 0 0 1]; c = [1 2]';

% Script-based tests: each section below is a separate test, and only the
% variables defined above are available in all sections

%% es_simulate_perturbation: no perturbation, no effect
for expansion = {'logarithmic', 'non-logarithmic'}
  for order = 1:2
    r = es_simulate_perturbation(N, W, ind_ext, ones(nr,1), s_ratio, result, expansion{1}, order);
    assert(max(abs(r.delta_c)) < 1e-10 && max(abs(r.delta_v)) < 1e-10);
  end
end

%% es_simulate_perturbation: small perturbation, all expansions agree
assert(max(abs(r1.delta_v)) > 1e-4);
assert(max(abs(r1.delta_c - r3.delta_c)) < 1e-3 * max(abs(r1.delta_c)) + 1e-4);
assert(max(abs(r2.delta_v - r4.delta_v)) < 1e-3);

%% es_simulate_perturbation: non-logarithmic first order is the linear response
delta_p = [result.u .* (u_ratio - 1); zeros(length(ind_ext),1)];
assert(max(abs(r3.delta_v - result.control.RJp_un * delta_p)) < 1e-10);

%% es_simulate_perturbation: stationary fluxes (logarithmic expansion)
assert(max(abs(N(setdiff(1:nm, ind_ext),:) * r2.v_new)) < 1e-10);

%% es_discover_perturbation: runs with known and unknown values
c_ratio = r1.c_ratio;
v_ratio = r1.v_ratio;
v_ratio(2) = nan;
[u_pred, s_pred, effect] = es_discover_perturbation(N, W, ind_ext, c_ratio, v_ratio, result, 'logarithmic', 1);
assert(isequal(size(u_pred), [nr 1]) && isequal(size(s_pred), [length(ind_ext) 1]));
assert(all(isfinite(effect.v_ratio)));

%% solve_quadratic: constraint is satisfied
x = solve_quadratic(A, c, B);
assert(max(abs(B * x - c)) < 1e-10);

%% solve_quadratic2: constraint is satisfied
D = zeros(2,3,3); D(1,:,:) = eye(3); D(2,:,:) = eye(3);
x = solve_quadratic2(A, c, B, D);
cpred = B * x + 1/2 * tensor_product(tensor_product(D, x, 2, 1), x, 2, 1);
assert(max(abs(cpred - c)) < 1e-6);
