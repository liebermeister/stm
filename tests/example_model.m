function [N, W, ind_ext, v, es_options, es_constraints] = example_model()

% [N, W, ind_ext, v, es_options, es_constraints] = example_model()
%
% Small network used in the tests:
%
%  S1             S5
%     \         /
%       S3 - S4
%     /         \
%  S2             S6
%
% S1, S2, S5, and S6 are external metabolites; S4 inhibits reaction 1.
% v is a stationary flux distribution; es_options and es_constraints are
% set up for sampling a model with these fluxes (half-saturated enzymes)

N = [ -1  0  0  0  0; ...
       0 -1  0  0  0; ...
       1  1 -1  0  0; ...
       0  0  1 -1 -1; ...
       0  0  0  1  0; ...
       0  0  0  0  1];

W = zeros(5,6);
W(1,4) = -1;

ind_ext = [1 2 5 6]';
v       = [1 1 2 1 1]';

[es_options, es_constraints] = es_default_options(N);
es_options.seed              = 1;
es_options.flag_second_order = 0;
es_constraints.v_fix         = v;
