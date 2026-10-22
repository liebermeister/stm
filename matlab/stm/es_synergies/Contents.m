% Functions for prediction of enzyme synergisms by elasticity sampling
%
% find_hyperantagonistic:
%     Check set of fitness 1st and 2nd derivatives
%     for existence of hyperantagonistic responses to gene decrease
%     (uses function check_hyperantagonistic)
%
% influence_anova:
%   Anova of influence (1st order) effects
%
% interaction_anova:
%   Anova of interaction (2nd order) effects
%
% interaction_statistics:
%   Run elasticity sampling repeatedly and compute statistics over influence and interaction values
%
% prism_clustering, sort_clusters, randomize_synergies, synergy_count_cycles:
%   Clustering of synergy networks and analysis of cycles
%
% For plotting interactions on a network, see synergy_network_plot (Metabolic Network Toolbox)
