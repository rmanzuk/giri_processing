input_folder = '/Users/rmanzuk/Desktop/bryo_project/classified_image_tifs';
output_folder = '/Users/rmanzuk/Desktop/bryo_project/binary_image_tifs';
filePattern = fullfile(input_folder, '*.tif');
tifFiles = dir(filePattern);
baseFileNames = natsortfiles({tifFiles.name});

for k = 1:numel(baseFileNames);
  fullFileName = fullfile(input_folder, baseFileNames{k});
  fprintf('Now reading %s\n', fullFileName);
  bw_image = imread(fullFileName);
  bw_image(bw_image == 1) = 0;
  bw_image(bw_image == 4) = 1;
  bw_image(bw_image == 3) = 1;
  bw_image(bw_image == 2) = 1;
  
  output_file_name = fullfile(output_folder, baseFileNames{k});
  imwrite(bw_image,output_file_name);
end