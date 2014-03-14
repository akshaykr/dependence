% Classifies images via distribution clustering

addpath ../common/
addpath ~/libs/libsvm/matlab/
clear all;
close all;

% filename = 'c12348.mat';
% filename = 'c128.mat';
% filename = 'c1to8.mat';
% filename = 'c1to8_n200_d5.mat';
% filename = 'C1to8_tr50te50_d10.mat';
filename = 'C1to8_tr50te50_d10_phow.mat';
load_data_from_file = true;
% load_data_from_file = false;

% Define the following constants
NUM_SIFT_DIMS = 128;
NUM_PCA_DIMS = 10;

% Parameters for the SVM Classifier
cv_cand_bws = [];
cv_cand_Cs = [0.001];
% svm_bandwidth = 0.01;
% svm_softmargin_penalty = 0.001;

% Prelims
if ~load_data_from_file
  num_train_images_per_cluster = 50;
  num_test_images_per_cluster = 50;

  categories = {'apples', 'cars', 'cows', 'cups', 'dogs', 'pears', ...
                'tomatoes', 'horses'};
%   categories = {'apples', 'cars', 'cows', 'cups', 'horses'};
%   categories = {'apples', 'cars', 'horses'};
%   categories = {'apples', 'horses'};

  % Obtain the Data
  num_categories = numel(categories);
  prepData;
  num_train_data = size(train_data, 1);
  num_test_data = size(test_data, 1);

  % Compute distance matrix
  all_data = [train_data; test_data];
  train_indices = 1:num_train_data;
  test_indices = (num_train_data + 1):(num_train_data + num_test_data);
  distance_matrix = computeL2DistanceMatrix(all_data);

  save(filename, 'num_train_images_per_cluster', ...
    'num_test_images_per_cluster', 'categories', 'num_categories', ...
    'num_train_data', 'num_test_data', 'train_indices', 'test_indices', ...
    'train_labels', 'test_labels', ...
    'distance_matrix');
else
  load(filename);
end

% Obtain the train-train and test-train distance matrices.
tr_distance_matrix = distance_matrix(1:num_train_data, 1:num_train_data);
tetr_distance_matrix = distance_matrix(num_train_data+1:end, 1:num_train_data); 

% Run the Classifier
model = distributionClassification(tr_distance_matrix, train_labels, ...
  cv_cand_bws, cv_cand_Cs);
% And make predictions
pred_labels = distributionClassificationPredict(model, tetr_distance_matrix);
tr_pred_labels = distributionClassificationPredict(model, tr_distance_matrix);

% Compute training and testing accuracy
tr_confmat_dg = confusion_matrix(tr_pred_labels, train_labels);
train_acc = accuracy_score(tr_pred_labels, train_labels);
fprintf('Train Results: confmat = %f, acc = %f\n', tr_confmat_dg, train_acc);
[conf_mat_diag, conf_mat] = confusion_matrix(pred_labels, test_labels);
test_acc = accuracy_score(pred_labels, test_labels);
fprintf('Test Results: confmat = %f, acc = %f\n', conf_mat_diag, test_acc);
conf_mat,

