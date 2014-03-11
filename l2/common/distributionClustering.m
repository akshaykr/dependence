function [labels, A] = distributionClustering(data, num_clusters, bandwidth)
% data is a cell array containing num_dists matrices. Each matrix of size
% nixnum_dims gives the points sampled from each distribution.
% bandwidth is the bandwidth for Spectral clustering.
% The kernel we use is K(p,q) = exp(-L2(p,q)^2/2h^2) where L2(p,q) is the L2
% divergence between p and q.
% Returns the labels and the similarity matrix.

  % Prelims
  num_dists = numel(data);
  
  % First compute the distance Matrix
  D = computeL2DistanceMatrix(data);
  A = exp(-D.^2/ (2*bandwidth^2));

  % Finally perform Spectral Clustering
  labels = spectralCluster(A, num_clusters);
end
