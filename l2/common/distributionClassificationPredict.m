function [pred_labels] = distributionClassificationPredict(model, ...
  test_train_dist_matrix)
% Inputs
% test_train_dist_matrix is a num_test_data x num_train_data matrix with the
% (i,j)^th element being the l2 distance between the ith test distribution and
% jth training distribution.
% Outputs
% model has 3 parameters. model.libsvm_model is the model obtained via svmtrain.
% model.bandwidth is the bandwidth for the SVM

  % Num test data
  num_test_data = size(test_train_dist_matrix, 1);

  % First construct the test-train kernel matrix
  tetr_kernel_matrix = exp(-test_train_dist_matrix.^2 / (2*model.bandwidth^2));
  k_tetr = [ones(num_test_data, 1) tetr_kernel_matrix];

  % Obtain labels
  pred_labels = svmpredict(ones(num_test_data, 1), ...
    k_tetr, model.libsvm_model, '-q');

end
