% Classifies images via distribution clustering

addpath ../common/
addpath ~/libs/libsvm/matlab/
addpath ../npdiv_export/
addpath ../npdiv_export/lib/
% clear all;
close all;

% Determine which tests to perform
performDistClass = true;
performBOWClass = true;
performBarnClass = true;

% filename = 'C1to8_tr50te50_d10_phow.mat'; % the old l2 estimator
% filename = 'C1to8_tr50te50_d10_phow_v2.mat'; % new l2 est with CV for bw
% filename = 'C1to8_tr50te50_d10_phow_sil.mat'; % new l2 est with silverman bw
filename = 'temp.mat';
% load_data_from_file = true;
load_data_from_file = false;

% Define the following constants
NUM_SIFT_DIMS = 384;
NUM_PCA_DIMS = 10;

% Parameters for the SVM Classifier
cv_cand_bws = [];
cv_cand_Cs = [];
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
  distance_matrix = computeL2DistanceMatrix(all_data, 1e-2);

  % Compute Barnabas's Stuff
  barn_Ds = sqrt( NPDivs(all_data, [], {@NPDivL22_RhoNu}, struct('k', '3')) );

  save(filename, 'num_train_images_per_cluster', ...
    'num_test_images_per_cluster', 'categories', 'num_categories', ...
    'num_train_data', 'num_test_data', 'train_indices', 'test_indices', ...
    'train_labels', 'test_labels', ...
    'distance_matrix', 'barn_Ds');
else
  load(filename);
end

if performDistClass
  model = distClassCVBW(train_data, train_labels, [], []);
  % Obtain the train-train and test-train distance matrices.
  tr_distance_matrix = distance_matrix(1:num_train_data, 1:num_train_data);
  tetr_distance_matrix = distance_matrix(num_train_data+1:end,1:num_train_data);

  % Run the Classifier
  model = distributionClassification(tr_distance_matrix, train_labels, ...
    cv_cand_bws, cv_cand_Cs);
  % And make predictions
  pred_labels = distributionClassificationPredict(model, tetr_distance_matrix);
  tr_pred_labels = distributionClassificationPredict(model, tr_distance_matrix);

  % Compute training and testing accuracy
  [tr_confmat_dg, tr_conf_mat] = confusion_matrix(tr_pred_labels, train_labels);
  train_acc = accuracy_score(tr_pred_labels, train_labels);
  fprintf('Train Results: confmat = %f, acc = %f\n', tr_confmat_dg, train_acc);
  [conf_mat_diag, conf_mat] = confusion_matrix(pred_labels, test_labels);
  test_acc = accuracy_score(pred_labels, test_labels);
  fprintf('Test Results: confmat = %f, acc = %f\n', conf_mat_diag, test_acc);
  conf_mat,
end

% Compare results with other methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. BOW
%%%%%%%%
if performBOWClass
  bow_class_model = bowClassification(train_data, train_labels, 300, []);
  [bow_train_preds] = bowClassificationPredict(bow_class_model, train_data);
  [bow_test_preds] = bowClassificationPredict(bow_class_model, test_data);
  % Compute training and testing accuracy
  [bow_tr_confmat_dg, bow_tr_conf_mat] = ...
    confusion_matrix(bow_train_preds, train_labels);
  bow_train_acc = accuracy_score(bow_train_preds, train_labels);
  fprintf('BOW Train Res: confmat = %f, acc = %f\n', ...
    bow_tr_confmat_dg, bow_train_acc);
  [bow_confmat_dg, bow_conf_mat] = confusion_matrix(bow_test_preds,test_labels);
  bow_test_acc = accuracy_score(bow_test_preds, test_labels);
  fprintf('Test Res: confmat = %f, acc = %f\n', bow_confmat_dg, bow_test_acc);
end

% 2. Barnabas & Co
%%%%%%%%%%%%%%%%%%
if performBarnClass
  % obtain the train-train and test-train distance matrices
%   Ds = sqrt( NPDivs(all_data, [], {@NPDivL22_RhoNu}, struct('k', '3')) );
  barn_tr_dist_matrix = barn_Ds(1:num_train_data, 1:num_train_data);
  barn_tetr_dist_matrix = barn_Ds(num_train_data+1:end, 1:num_train_data);

  % Run classifier
  barnModel = distributionClassification(barn_tr_dist_matrix, ...
    train_labels, cv_cand_bws, cv_cand_Cs);
  % And make predictions
  barn_pred_labels = distributionClassificationPredict(barnModel, ...
                      barn_tetr_dist_matrix);
  barn_tr_pred_labels = distributionClassificationPredict(barnModel, ...
                          barn_tr_dist_matrix);

  % Compute training and testing accuracy
  [barn_tr_confmat_diag, barn_tr_conf_mat] = ...
    confusion_matrix(barn_tr_pred_labels, train_labels);
  barn_train_acc = accuracy_score( barn_tr_pred_labels, train_labels);
  fprintf('Barn Train Res: confmat = %d, acc = %f\n', ...
    barn_tr_confmat_diag, barn_train_acc);
  % testing
  [barn_te_confmat_diag, barn_te_conf_mat] = ...
    confusion_matrix(barn_pred_labels, test_labels);
  barn_test_acc = accuracy_score( barn_pred_labels, test_labels);
  fprintf('Barn Test Res: confmat = %d, acc = %f\n', ...
    barn_te_confmat_diag, barn_test_acc);
end

