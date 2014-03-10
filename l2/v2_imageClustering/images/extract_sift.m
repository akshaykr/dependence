% Extracts Sift features from all the images in the current directory.
% Uses the VLFeat library.
% Specify the filename and the directory as a workspace variable when running.

categories = {'apples', ...
              'cars', ...
              'cows', ...
              'cups', ...
              'dogs', ...
              'pears', ...
              'tomatoes' ...
              'horses'};

for cat_iter = 1:numel(categories)
  dirname = categories{cat_iter};
  dir_search_name = sprintf('%s/*.png', dirname);

  imfiles = dir(dir_search_name);
  num_images = numel(imfiles);
  images = cell(num_images, 1);

  for i = 1:num_images
    imfile_name = sprintf('%s/%s', dirname, imfiles(i).name);
    img = imread(imfile_name);
    [~, features] = vl_sift(single(rgb2gray(img)));
    images{i} = double(features');
  end

  % Save the file
  save_file_name = sprintf('%s.mat', dirname);
  save(save_file_name, 'images');
end

