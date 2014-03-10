% Clusters the images in the eth-80 dataset

addpath ../common/
clear all;
close all;

% Define the following constants
NUM_SIFT_DIMS = 128;
NUM_PCA_DIMS = 4;
SC_BANDWIDTH = 0.001; % bandwidth for spectral clustering

num_images_per_cluster = 60;

categories = {'apples', ...
%               'cars', ...
%               'cows', ...
%               'cups', ...
%               'dogs', ...
%               'pears', ...
%               'tomatoes' ...
              'horses'};

% Obtain the Data
num_clusters = numel(categories);
prepData;

% Cluster and report results
fprintf('Clustering the Data\n');
[pred_labels, A] = distributionClustering(data, num_clusters, SC_BANDWIDTH);

% Compute accuracy
fprintf('Avg_diag of confusion matrix: %f\n', ...
  confusion_matrix(pred_labels, true_labels));

