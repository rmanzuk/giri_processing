input_folder = '/Users/rmanzuk/Desktop/bryo_project/binary_image_tifs';
output_folder = '/Users/rmanzuk/Desktop/bryo_project/binary_image_tifs_zero_filled';
filePattern = fullfile(input_folder, '*.tif');
tifFiles = dir(filePattern);
baseFileNames = natsortfiles({tifFiles.name});

for k = 1:numel(baseFileNames)
  fullFileName = fullfile(input_folder, baseFileNames{k});
  fprintf('Now reading %s\n', fullFileName);
  bw_image = imread(fullFileName);
  imshow(double(bw_image));
  [col_coord,row_coord] = ginput(2);
  bw_image(1:row_coord(1),:) = 0;
  bw_image(row_coord(2):end,:) = 0;
  bw_image(:,1:col_coord(1)) = 0;
  bw_image(:,col_coord(2):end) = 0;
  

  output_file_name = fullfile(output_folder, baseFileNames{k});
  imwrite(bw_image,output_file_name);
end