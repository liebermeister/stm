function [N, W, ind_ext, v, stm_options, stm_constraints] = example_model()

% [N, W, ind_ext, v, stm_options, stm_constraints] = example_model()
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
% v is a stationary flux distribution; stm_options and stm_constraints are
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

[stm_options, stm_constraints] = stm_default_options(N);
stm_options.seed              = 1;
stm_options.flag_second_order = 0;
stm_constraints.v_fix         = v;
