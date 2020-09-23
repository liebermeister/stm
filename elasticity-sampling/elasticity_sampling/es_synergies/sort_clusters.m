function [new_cluster_indices,new_cluster_interactions] = sort_clusters(network,cluster_indices,cluster_interactions)

[nm,nr] = network_numbers(network);

nc = max(cluster_indices);

for it = 1:nc,
 ypos(it) = mean(network.graphics_par.x(2,nm+find(cluster_indices==it)));
end

[yposnew,order] = sort(ypos);

for it = 1:nc,
 new_cluster_indices(find(cluster_indices==order(it))) = it;
end

if nargout >1,
  new_cluster_interactions = cluster_interactions(order,order);
end