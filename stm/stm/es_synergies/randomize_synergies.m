function E_rnd = randomize_synergies(E)

% E_rnd = randomize_synergies(E)
%
% Randomisation of the interaction values in the upper triangle (without diagonal)
% of a symmetric matrix E; the result is symmetric and has a zero diagonal

n   = size(E,1);
ind = find(triu(ones(n),1));
ind_perm = ind(randperm(length(ind)));
E_rnd = zeros(n);
E_rnd(ind_perm) = E(ind);
E_rnd = E_rnd + E_rnd';
