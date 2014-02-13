% Distribution clustering for 2 clusters

% Constants
NUM_DIMS = 10;
NUM_PTS = 100; % The number of points to draw from from each distribution
NUM_DISTS = 100; % The number of distributions
% Make sure that NUM_PTS * NUM_DISTS * NUM_DIMS < 50,000.

genData;

% Parameters for clustering 
bandwidth = 0.05;

% Now construct Similarity matrix for each distribution
A = eye(NUM_DISTS);
for i = 1:NUM_DISTS
  fprintf('Computing similarities for i = %d\n', i);
  for j = (i+1):NUM_DISTS
    l2ij = l2Divergence(Data{i}, Data{j});
    k2ij = exp(- l2ij^2/ (2*bandwidth^2) );
    A(i,j) = k2ij;
    A(j,i) = k2ij;
  end
end

% Peform Spectral Clustering
labels = spectralCluster(A, 2);

% Compare clusering
fprintf('Avg_diag of confusion matrix: %f\n', ...
  confusion_matrix(labels, mixture_idxs));
% [labels, mixture_idxs],

