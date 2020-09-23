function E_rnd = randomize_interaction(E)

% E_rnd = randomize_interaction(E)
% randomisation of interaction values with the upper triangle

n   = size(E,1);
ind = find(triu(ones(n)));
ind_perm = ind(randperm(length(ind)));
E_rnd = zeros(n);
E_rnd(ind_perm) = E(ind);
E_rnd = E_rnd + E_rnd';
