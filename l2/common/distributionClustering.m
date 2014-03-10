function [labels, A] = distributionClustering(data, num_clusters, bandwidth)
% data is a cell array containing num_dists matrices. Each matrix of size
% nixnum_dims gives the points sampled from each distribution.
% bandwidth is the bandwidth for Spectral clustering.
% The kernel we use is K(p,q) = exp(-L2(p,q)^2/2h^2) where L2(p,q) is the L2
% divergence between p and q.
% Returns the labels and the similarity matrix.

  % Prelims
  num_dists = numel(data);
  
  % First compute the similarity matrix
  A = zeros(num_dists);
  for i = 1:num_dists
    if (num_dists > 40) && (mod(i,10) == 0)
      % then report progress
      fprintf('Computing L2 between %dth distribution and others\n', i);
    end
    for j = (i+1):num_dists
      l2ij = l2Divergence(data{i}, data{j});
      k2ij = exp(-l2ij^2/ (2*bandwidth^2) );
      A(i,j) = k2ij;
      A(j,i) = k2ij;
    end
  end

  % Finally perform Spectral Clustering
  labels = spectralCluster(A, num_clusters);
end
