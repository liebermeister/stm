function n_min = fdr_check_permutations(p, n_per, fdr)

% FDR_CHECK_PERMUTATIONS - Warn if there are too few permutations for a false discovery rate
%
% n_min = fdr_check_permutations(p, n_per, fdr)
%
% With n_per permutations, the smallest two-sided p value is 2/(n_per+2) (see pvalue_from_countnumber
% and fdr_two_sided). In the Benjamini-Hochberg procedure with n tests, such p values
% are significant only if there are at least n_min = n * 2/(n_per+2) / fdr of them.
%
% p:     p values for all tests (nan values are not counted as tests)
% n_per: number of permutations
% fdr:   false discovery rate
%
% n_min: smallest number of significant results that can be detected.
%        If n_min > 1, a warning is shown.

n_tests = sum(isfinite(p(:)));
n_min   = max(1, ceil(n_tests * 2 / (n_per + 2) / fdr));

if n_min > 1
  warning('%d permutations are too few for %d tests with FDR %g: fewer than %d significant results cannot be detected. Use at least %d permutations.', ...
          n_per, n_tests, fdr, n_min, ceil(2 * n_tests / fdr) - 2);
end
