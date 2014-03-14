% Extracts Sift features from all the images in the current directory.
% Uses the VLFeat library.
% Specify the filename and the directory as a workspace variable when running.

addpath ~/libs/vlfeat/toolbox/
vl_setup;

categories = {'apples', ...
              'cars', ...
              'cows', ...
              'cups', ...
              'dogs', ...
              'pears', ...
              'tomatoes', ...
              'horses'};

for cat_iter = 1:numel(categories)
  dirname = categories{cat_iter};
  dir_search_name = sprintf('%s/*.png', dirname);

  imfiles = dir(dir_search_name);
  num_images = numel(imfiles);
  images = cell(num_images, 1);

  for i = 1:num_images
%     fprintf('Processing %d/%d\n', i, num_images);
    imfile_name = sprintf('%s/%s', dirname, imfiles(i).name);
    img = imread(imfile_name);
    [~, features] = vl_phow(single(img), 'ContrastThreshold', 0, ...
      'color', 'rgb', 'FloatDescriptors', true, 'step', 6, 'sizes', [6 9 12]);
%     [~, features] = vl_dsift(single(rgb2gray(img)), 'step', 4, 'size', 5);
    images{i} = double(features');
  end

  % Save the file
  save_file_name = sprintf('%s.mat', dirname);
  save(save_file_name, 'images');
end

