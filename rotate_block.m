myFolder = '/Users/rmanzuk/Desktop/nevada_jan_2019/sm_117_71_reconstruction/cleaned_binary_tifs';
filePattern = fullfile(myFolder, '*.tiff');
tifFiles = dir(filePattern);
output_folder = '/Users/rmanzuk/Desktop/nevada_jan_2019/sm_117_71_reconstruction/1st_rotation_binaries';

baseFileName = tifFiles(1).name;
fullFileName = fullfile(myFolder, baseFileName);
image = imread(fullFileName);
bw_image = im2bw(image,0.5);
imshow(bw_image);
[crop_spot_x,crop_spot_y] = ginput;


for k = 1:length(tifFiles)
  baseFileName = tifFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf('Now reading %s\n', fullFileName);
  image = imread(fullFileName);
  bw_image = im2bw(image,0.5);
  crop_image = bw_image(:,round(crop_spot_x):end);
  rot_im = imrotate(crop_image,307,'loose');
  fullFileName2 = fullfile(output_folder, baseFileName);
  imwrite(rot_im, fullFileName2);
end

   