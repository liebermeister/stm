% Simulation of perturbed metabolic systems based on elasticity sampling:
% perturbations by differential expression and varying external metabolite levels
% 
% Functions
%  es_simulate_perturbation    Predict the effects of perturbations
%  es_discover_perturbation    Estimate most likely perturbation behind a given effect 
%
% Same functions, embedded in a loop for sampling:
%  es_simulate_perturbation_sample
%  es_discover_perturbation_sample