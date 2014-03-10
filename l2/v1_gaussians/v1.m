% Distribution clustering for 2 clusters
addpath ../common/

% Constants
NUM_DIMS = 2;
NUM_PTS = 80; % The number of points to draw from from each distribution
NUM_DISTS = 50; % The number of distributions
% Make sure that NUM_PTS * NUM_DISTS * NUM_DIMS < 50,000.

genData;

% Parameters for clustering 
bandwidth = 1;

% Peform Spectral Clustering
fprintf('Performing Distribution Clustering. Will take some time ... \n');
labels = distributionClustering(Data, 2, bandwidth);

% Compare clusering
fprintf('Avg_diag of confusion matrix: %f\n', ...
  confusion_matrix(labels, mixture_idxs));

