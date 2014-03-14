function model = bowClassification(trainDistros, Ytr, vocabSize)
% I am using vl_feat here so that needs to be initialized when running this code

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
  model.codeBook = kmeans(all_features, vocabSize);
  % Create a function handle to retrieve the cluster index
  model.getCodeBookIdx = @(arg) retrieveKMeansClusterIds(arg, model.code_book);
  % and code book representation of a set of distributions
  model.obtainCodeBookRepr = @(arg) obtainCodeBookRepr(arg, model.codeBook);
  
  % Obtain code book representation of the training data
  Xtr = model.obtainCodeBookRepr(trainDistros);
  % normalize the data to have zero mean and unit variance
  cb_means = mean(Xtr)';
  cb_stds = std(Xtr)'; cbs_stds = cb_stds + double(cbs_stds == 1);
  nXtr = bsxfun(@rdivide, bsxfun(@minus, Xtr, cb_means'), cb_stds');
  % Store these results in model
  model.cb_means = cb_means;
  model.cb_stds = cb_stds;
  model.obtainNormCodeBookRepr = @(arg) bsxfun(@rdivide, ...
    bsxfun(@minus, obtainCodeBookRepr(arg, model.codeBook), cb_means') ...
    cb_stds');

  % Now finally, train the svm


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
    X(i,:) = retrieveKMeansClusterIds(distros{i}, code_book)';
  end
end

