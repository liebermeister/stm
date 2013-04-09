function p = pvalue_from_countnumber(value,value_list)

% formula stems from the mean of a beta distribution
% posterior of p based on binominal distribution for
% value<=value_list and flat prior for p

N = length(value_list);
n = length(find(value<=value_list));
p = (n+1)/(N+2);