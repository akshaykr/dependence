function model = bowClassification(trainDistros, Ytr, vocabSize, ...
  cv_cand_Cs)
% I am using vl_feat here so that needs to be initialized when running this code

  DFLT_NUM_CANDIDATES = 15; % number of candidates for the bandwidth and C 
  NUM_KFOLDCV_PARTITIONS = 10;

  % First accumulate all features
  num_feat_dims = size(trainDistros{1}, 2);
  num_train_dists = numel(trainDistros);
  class_idxs = unique(Ytr);
  num_classes = size(Ytr, 2);
  all_features = zeros(0, num_feat_dims);
  total_num_features = 0;
  for train_idx = 1:num_train_dists
    total_num_features = total_num_features + size(trainDistros{train_idx}, 1);
    all_features = [all_features; trainDistros{train_idx}];
  end

  % Now obtain the centres
  [~, codeBook] = kmeans(all_features, vocabSize);
  model.codeBook = codeBook;
  % Create a function handle to retrieve the cluster index
  model.getCodeBookIdx = @(arg) retrieveKMeansClusterIds(arg, model.codeBook);
  % and code book representation of a set of distributions
  model.obtainCodeBookRepr = @(arg) obtainCodeBookRepr(arg, model.codeBook);
  
  % Obtain code book representation of the training data
  Xtr = model.obtainCodeBookRepr(trainDistros);
  normalize = true;
%   normalize = false;
  if normalize
  % normalize the data to have zero mean and unit variance
    cb_means = mean(Xtr)';
    cb_stds = std(Xtr)'; cb_stds = cb_stds + double(cb_stds == 0);
    nXtr = bsxfun(@rdivide, bsxfun(@minus, Xtr, cb_means'), cb_stds');
    % Store these results in model
    model.cb_means = cb_means;
    model.cb_stds = cb_stds;
    model.obtainNormCodeBookRepr = @(arg) bsxfun(@rdivide, ...
      bsxfun(@minus, model.obtainCodeBookRepr(arg), cb_means'), ...
      cb_stds');
  else
    model.obtainNormCodeBookRepr = model.obtainCodeBookRepr;
    nXtr = Xtr;
  end

  % Now  perform CV over the soft margin constants
  % shuffle the data
  shuffle_order = randperm(num_train_dists);
  shuffled_nXtr = nXtr(shuffle_order, :);
  shuffled_Ytr = Ytr(shuffle_order, :);
  % determine the Soft margin candidates
  if isempty(cv_cand_Cs)
    cv_cand_Cs = logspace(-3, 1, DFLT_NUM_CANDIDATES)';
  end
  % Perform CV 
  best_cv_acc = 0;
  for c_iter = 1:numel(cv_cand_Cs)
    libsvm_args = sprintf('-c %f -v %d -q', cv_cand_Cs(c_iter), ...
      NUM_KFOLDCV_PARTITIONS);
    curr_cv_acc = svmtrain(shuffled_Ytr, shuffled_nXtr, libsvm_args);
    if curr_cv_acc >= best_cv_acc
      best_cv_acc = curr_cv_acc;
      best_C = cv_cand_Cs(c_iter);
    end
  end
  % Print out values picked
  fprintf('Picked C = %0.4f, cv_acc = %0.4f\n', best_C, best_cv_acc);

  % Finally train the whole thing again
  libsvm_args = sprintf('-c %f -q', best_C);
  libsvm_model = svmtrain(Ytr, nXtr, libsvm_args);
  libsvm_model,
  % test on the data and report it
  [train_preds, acc] = svmpredict(Ytr, nXtr, libsvm_model, '-q'),
%   accuracy_score(train_preds, Ytr),
  model.libsvm_model = libsvm_model;
  model.C = best_C;

end

% This function retrieves the cluster index of a point in K-means
function [clust_idxs] = retrieveKMeansClusterIds(X, centres)
% X is an nxd data matrix. centres is a num_centresxd data matrix
  [~, clust_idxs] = min(dist2(X, centres), [], 2);
end

% This function returns the code book representation of a cell of distributions.
function X = obtainCodeBookRepr(distros, code_book)
  num_distros = numel(distros);
  vocab_size = size(code_book, 1);
  X = zeros(num_distros, vocab_size);
  for i = 1:num_distros
    codeBookIds = retrieveKMeansClusterIds(distros{i}, code_book);
    t = tabulate(codeBookIds);
    X(i,:) = [t(:,2); zeros(vocab_size - size(t, 1), 1)]';
  end
end

