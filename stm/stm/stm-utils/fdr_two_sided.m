function significant = fdr_two_sided(p, fdr)

% FDR_TWO_SIDED - Two-sided significance from one-sided p values
%
% significant = fdr_two_sided(p, fdr)
%
% p:   one-sided p values (small p: value is high, large p: value is low)
% fdr: false discovery rate for all two-sided tests together
%
% significant: 1 (significantly high), -1 (significantly low), 0 (not significant or p = nan)

significant = double(multiple_testing_fdr(2 * min(p,1-p), fdr)) .* sign(0.5 - p);
significant(isnan(p)) = 0;
