% order=sort_by_clustering(X)
% 
% sort columns of matrix X

function [order,h,l] = sort_by_clustering(X)

%l=linkage(pdist(X','cityblock'));
l=linkage(pdist(X','correlation'));
figure(10000);h = dendrogram(l,0);
order = str2num(get(gca,'XTickLabel'));
