function[n_cycles,frac_pos,mean_n_cycles,std_n_cycles,pvalue_n_cycles,pvalue_frac_cycles] = interaction_count_cycles(E,n_rand)

% E: epistasis matrix;

E_sign_triu    = sign(triu(E,1));
E_nonzero_triu = double(E_sign_triu~=0);

% n_cycles: for each node n (row) and each cycle length l
% number of cycles of length l starting and ending at n
% edges may only be used once in a cycle

% d_cycles: for each node n (row) and each cycle length l
% similar: number of positive minus number of negative cycles
% the sign of a cycle is the product of the signs of its edges

n_cycles=[ 0 0];
d_cycles=[ 0 0];

for it = 3:10,
  n_cycles(it) = sum(sum(E_nonzero_triu .* E_nonzero_triu^(it-1) ));
  d_cycles(it) = sum(sum(E_sign_triu    .* E_sign_triu^(it-1)    ));
end

pos = 0.5 * (n_cycles + d_cycles);
neg_cycles = 0.5 * (n_cycles - d_cycles);

frac_pos = pos./(n_cycles+10^-10);
frac_pos(n_cycles==0) = nan;

if exist('n_rand','var'),
  
  for it = 1:n_rand,
    E_rnd = randomize_interaction(E);
    [n_cycles_rnd(it,:),frac_pos_rnd(it,:)] = interaction_count_cycles(E_rnd);
  end
  mean_n_cycles      = nanmean(n_cycles_rnd);
  std_n_cycles       = nanstd(n_cycles_rnd);
  pvalue_n_cycles    = nanmean(repmat(n_cycles,n_rand,1)>n_cycles_rnd);
  pvalue_frac_cycles = nanmean(repmat(frac_pos,n_rand,1)>frac_pos_rnd);
end
