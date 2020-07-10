myFolder = '/Users/rmanzuk/Desktop/nevada_jan_2019/sm_117_71_reconstruction/bedding_plane_pass2'
filePattern = fullfile(myFolder, '*.tiff');
tifFiles = dir(filePattern);
baseFileNames = natsortfiles({tifFiles.name});
num_components = zeros(length(tifFiles),1);
component_mean = zeros(length(tifFiles),1);
depth = [1:length(tifFiles)];
depth = (depth .* 20)./10000;
for k = 1:numel(baseFileNames);
  fullFileName = fullfile(myFolder, baseFileNames{k});
  fprintf('Now reading %s\n', fullFileName);
  bw_image = imread(fullFileName);
  bw_image2 = bwareafilt(bw_image,[15000,80000000]);
  CC = bwconncomp(bw_image2);
  numPixels = cellfun(@numel,CC.PixelIdxList);
  component_mean(k) = mean(numPixels);
  mean_comp_cm = (component_mean .* 400) ./ 100000000;
  num_components(k) = CC.NumObjects;
  plot(depth,num_components,'LineStyle','-');
  drawnow;
end

