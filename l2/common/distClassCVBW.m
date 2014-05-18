function [model] = distClassCVBW( ...
  train_data, train_labels, cv_cand_smoothing_bws, cv_cand_Cs)
% function [model] = distributionClassification( ...
%   train_distance_matrix, train_labels, cv_cand_bws, cv_cand_Cs)
% This is a wrapper for distribution Classification but also does the CV on
% the bandwidths.
% train_distance_matrix is the matrix of distances between each distribuiton.
% Returns model: a struct containing the libsvm_model, and the optimal params.

  DFLT_NUM_CANDIDATES = 10; % number of candidates for the bandwidth and C 
  NUM_KFOLDCV_PARTITIONS = 10;

  % prelims
  num_train_data = numel(train_data);

  % determine the bandwidth candidates
  if isempty(cv_cand_smoothing_bws)
    cv_cand_smoothing_bws = logspace(10^-3, 10, DFLT_NUM_CANDIDATES)';
  end

  % determine the soft margin constant candidates.
  if isempty(cv_cand_Cs)
    cv_cand_Cs = logspace(-3, 1, DFLT_NUM_CANDIDATES)';
  end

  num_sm_bw_cands = numel(cv_cand_smoothing_bws);
  num_C_cands = numel(cv_cand_Cs);

  % Shuffle the data for K-Fold cross validation
  data_ordering = randperm(num_train_data);
  shuffled_tr_data = train_data(data_ordering);
  shuffled_tr_labels = train_labels(data_ordering, :);

  % Use K-fold cross validation to pick k
  best_cv_acc = 0;
  for i = 1:num_sm_bw_cands % iterate through cv_cand_bws

    % First prepare the tr_tr_kernel matrix
%     shuffled_trtr_matrix = computeL2DistanceMatrix(shuffled_tr_data, ...
%       cv_cand_smoothing_bws(i) );
    R = rand(num_train_data); shuffled_trtr_matrix = R * R';

    % Prep arguments for svmtrain
    train_kernel_mat = exp(-shuffled_trtr_matrix.^2 / 2 );
    k_train = [(1:num_train_data)' train_kernel_mat];

    for j = 1:num_C_cands % iterate through cv_cand_Cs

      % -t 4: we are using a custom kernel. -c C is the soft margin constant
      % -v K: use K fold CV, -q: quiet mode
      libsvm_args = sprintf('-t 4 -c %f -v %d -q', ...
        cv_cand_Cs(j), NUM_KFOLDCV_PARTITIONS);
      % Run svmtrain
      curr_cv_acc = svmtrain(shuffled_tr_labels, k_train, libsvm_args);
%       fprintf('CV Results (h,C) = (%e, %f) : %f\n', ...
%         cv_cand_bws(i), cv_cand_Cs(j), curr_cv_acc);

      % Determine if this is the best value
      if best_cv_acc < curr_cv_acc;
        best_cv_acc = curr_cv_acc;
        best_bw = cv_cand_bws(i);
        best_C = cv_cand_Cs(j);
      end

    end
  end

  % Print out the values you picked
  fprintf('bandwidth candidates: %s\n', mat2str(cv_cand_smoothing_bws));
  fprintf('Soft-margin-const candidates: %s\n', mat2str(cv_cand_Cs));
  fprintf('Picked (bw,C) = (%.e, %.4f). CV accuracy: %.4f\n', ...
          best_bw, best_C, best_cv_acc);
    
  % Finally, retrain the SVM using all parameter and the data
  train_distance_matrix = computeL2DistanceMatrix(train_data, best_bw);
  train_kernel_mat = exp(-train_distance_matrix.^2 / 2);
  k_train = [(1:num_train_data)' train_kernel_mat];
  libsvm_args = sprintf('-t 4 -c %f -q', best_C);
  libsvm_model = svmtrain(train_labels, k_train, libsvm_args);

  % Return the structure
  model.bandwidth = best_bw;
  model.C = best_C;
  model.libsvm_model = libsvm_model;

end

