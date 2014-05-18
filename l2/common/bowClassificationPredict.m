function [preds] = bowClassificationPredict(model, test_distros)
  % obtain the normalized bow Representation
  Xte = model.obtainNormCodeBookRepr(test_distros);
  num_test_data = size(Xte, 1);
  % Make predictions using libsvm
  preds = svmpredict(ones(num_test_data, 1), Xte, model.libsvm_model, '-q');
end

