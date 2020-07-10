input_folder = '/Users/rmanzuk/Desktop/nevada_jan_2019/sm_117_71_reconstruction/bedding_plane_pass';
filePattern = fullfile(input_folder, '*.tiff');
tifFiles = dir(filePattern);
output_folder = '/Users/rmanzuk/Desktop/nevada_jan_2019/sm_117_71_reconstruction/bedding_plane_pass2';

for k = 3716:length(tifFiles);
  baseFileName = tifFiles(k).name;
  fullFileName = fullfile(input_folder, baseFileName);
  fprintf('Now reading %s\n', fullFileName);
  image = imread(fullFileName);

  bw_image = im2bw(image);
  bw_inv_image = imcomplement(bw_image);

  check_column = all(bw_inv_image);

  for i = 1:size(bw_inv_image,2)
      if check_column(i) == 1
          bw_inv_image(:,i) = 0;
      end
  end
  fullFileName2 = fullfile(output_folder, baseFileName);
  imwrite(bw_inv_image, fullFileName2);
  
end


        

