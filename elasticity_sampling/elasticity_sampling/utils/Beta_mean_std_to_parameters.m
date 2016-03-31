function [Beta_A_alpha, Beta_A_beta] = Beta_mean_std_to_parameters(x_mean,x_std)

%[Beta_A_alpha, Beta_A_beta] = Beta_mean_std_to_parameters(x_mean,x_std)

% Based on formulae for computing distributions parameters from mean (x_mean) and std_dev (x_std)
% "method of moments"; works for matrices

dum          = [ x_mean .* [1- x_mean] ] ./ [x_std.^2] -1;
Beta_A_alpha = x_mean .* dum;
Beta_A_beta  = [1-x_mean] .* dum;
