% This script will load the data from the images directory

% First load all images and obtain the principal components
total_num_images = 0;
all_vectors = zeros(0, NUM_SIFT_DIMS);
for i = 1:num_categories
  % load data
  fprintf('Loading %s ... ', categories{i});
  file_name = sprintf('images/%s.mat', categories{i});
  load(file_name);
  num_curr_images = numel(images);
  total_num_images = total_num_images + num_curr_images;
  fprintf('%d images, ', num_curr_images);
  % Add all data here to all_vectors
  num_vecs_in_categ = 0;
  for j = 1:num_curr_images %images is the name of the struct
    all_vectors = [all_vectors; images{j}];
    num_vecs_in_categ = num_vecs_in_categ + size(images{j}, 1);
  end
  fprintf('%d features\n', num_vecs_in_categ);
end
fprintf('Total # images: %d\n', total_num_images);

%% Normalize the data
data_mean = mean(all_vectors)';
data_std = std(all_vectors)';
% In case there are redundant features use std = 1;
data_std = data_std + double((data_std==0));
all_vectors = bsxfun(@minus, all_vectors, data_mean') / diag(data_std);
% Now use SVD. Plot the singular values
fprintf('Running PCA ...\n');
[~, S, V] = svd(all_vectors, 0); % use thin svd
figure; plot(diag(S).^2); title('Eigenvalues of the data matrix');
clear all_vectors;

% Use only the first NUM_PCA_DIMS principal components
% Load only num_train_images_per_cluster images per data
T = V(:, 1:NUM_PCA_DIMS);
% Reread the data and and compress them
train_labels = zeros(0, 1);
test_labels = zeros(0, 1);
train_data = cell(num_train_images_per_cluster * num_categories, 1);
test_data = cell(num_test_images_per_cluster * num_categories, 1);
train_data_counter = 0;
test_data_counter = 0;

for i = 1:num_categories
  % load data
  fprintf('Loading %s \n', categories{i});
  file_name = sprintf('images/%s.mat', categories{i});
  load(file_name);
  num_curr_images = numel(images);
  % Now select which images to choose randomly
  shuffled_images = randperm(num_curr_images);
  train_images = shuffled_images(1:num_train_images_per_cluster);
  test_images = shuffled_images(num_train_images_per_cluster+1: ...
    num_train_images_per_cluster+num_test_images_per_cluster);

  % Now add the train data after transforming them 
  for img_idx = train_images %images is the name of the struct
    train_data_counter = train_data_counter + 1;
    rescaled_img = bsxfun(@minus, images{img_idx}, data_mean') / diag(data_std);
    train_data{train_data_counter} = rescaled_img * T;
  end
  % ADd the test data
  for img_idx = test_images
    test_data_counter = test_data_counter + 1;
    rescaled_img = bsxfun(@minus, images{img_idx}, data_mean') / diag(data_std);
    test_data{test_data_counter} = rescaled_img * T;
  end

  % Add curr_categ_images to data
  % also save the labels.
  train_labels = [train_labels; i * ones(num_train_images_per_cluster, 1)];
  test_labels = [test_labels; i * ones(num_test_images_per_cluster, 1)];
end

