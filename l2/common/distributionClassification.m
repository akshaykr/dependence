function [model] = distributionClassification(...
  train_distance_matrix, train_labels, cv_cand_bws, cv_cand_Cs);
% train_distance_matrix is the matrix of distances between each distribuiton.
% Returns model: a struct containing the libsvm_model, and the optimal params

  DFLT_NUM_CANDIDATES = 15; % number of candidates for the bandwidth and C 
  NUM_KFOLDCV_PARTITIONS = 10;

  % prelims
  num_train_data = size(train_distance_matrix, 1);

  % determine the bandwidth candidates
  if isempty(cv_cand_bws)
    u = sort(unique(train_distance_matrix)); u = u(2:end);
    bw_min = 0.5 * u(1); 
    bw_max = 2 * u(end);
    cv_cand_bws = logspace(log10(bw_min), log10(bw_max), DFLT_NUM_CANDIDATES);
  end

  % determine the soft margin constant candidates.
  if isempty(cv_cand_Cs)
    cv_cand_Cs = logspace(-3, 1, DFLT_NUM_CANDIDATES);
  end

  num_bw_cands = numel(cv_cand_bws);
  num_C_cands = numel(cv_cand_Cs);

  % Use K-fold cross validation to pick k
  best_cv_acc = 0;
  for i = 1:num_bw_cands % iterate through cv_cand_bws
    for j = 1:num_C_cands % iterate through cv_cand_Cs

      % Prep arguments for svmtrain
      train_kernel_mat = exp(-train_distance_matrix.^2 / ...
        (2 * cv_cand_bws(i)^2) );
      % shuffle train_kernel_mat
      shuffle_order = randperm(num_train_data);
      train_kernel_mat = train_kernel_mat(shuffle_order, shuffle_order);

      k_train = [(1:num_train_data)' train_kernel_mat];
      % -t 4: we are using a custom kernel. -c C is the soft margin constant
      % -v K: use K fold CV, -q: quiet mode
      libsvm_args = sprintf('-t 4 -c %f -v %d -q', ...
        cv_cand_Cs(j), NUM_KFOLDCV_PARTITIONS);
      % Run svmtrain
      curr_cv_acc = svmtrain(train_labels, k_train, libsvm_args);

      % Determine if this is the best value
      if best_cv_acc < curr_cv_acc;
        best_cv_acc = curr_cv_acc;
        best_bw = cv_cand_bws(i);
        best_C = cv_cand_Cs(j);
      end

    end
  end

  % Print out the values you picked
  fprintf('bandwidth candidates: %s\n', mat2str(cv_cand_bws));
  fprintf('Soft-margin-const candidates: %s\n', mat2str(cv_cand_Cs));
  fprintf('Picked (bw,C) = (%.e, %.4f). CV accuracy: %.4f\n', ...
          best_bw, best_C, best_cv_acc);
    
  % Finally, retrain the SVM using all parameter and the data
  train_kernel_mat = exp(-train_distance_matrix.^2 / (2 * best_bw^2) );
  k_train = [(1:num_train_data)' train_kernel_mat];
  libsvm_args = sprintf('-t 4 -c %f -q', best_C);
  libsvm_model = svmtrain(train_labels, k_train, libsvm_args);

  % Return the structure
  model.bandwidth = best_bw;
  model.C = best_C;
  model.libsvm_model = libsvm_model;

end
