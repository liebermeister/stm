function significant = multiple_testing_fdr(all_p, false_discovery_rate)

% significant = multiple_testing_fdr(all_p, false_discovery_rate)
%
% Significance (yes/no) for multiple tests at a given false discovery rate
%
% Formula from http://en.wikipedia.org/wiki/False_discovery_rate, 
%         assuming independent tests!!!
%
% all_p (vector or matrix) contains the p values for all tests
%                          'nan' values are ignored

p           = all_p(isfinite(all_p)); 
if sum(p)==0, error('All p values are zero'); end
if length(p)==0, error('No finite p values found'); end
n_tests     = length(p);
alpha       = false_discovery_rate * 2 * n_tests/(n_tests+1);
sorted_p    = sort(p);
n_accept    = max(find( ...
     sorted_p <=  ( (1:n_tests)' / n_tests * alpha )...
    ) );
if length(n_accept),
  p_accept    = sorted_p(n_accept);
else
  p_accept    = 0;
end
significant = (all_p <= p_accept);

% plot(sorted_p); hold on; plot((1:n_tests)' / n_tests * alpha,'r'); hold off;
