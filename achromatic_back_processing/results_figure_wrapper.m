% Runs analyses on several mud stones from an Ardmore plate. Okay result,
% not used in the final paper.
%
% R. A. Manzuk 03/25/2021
%% load everything
im_dir = '/Users/ryan/Desktop/achromatic_project/p28_72mm';

file_pattern = fullfile(im_dir, '*.tif');
tifs = dir(file_pattern);
n_imgs = numel(tifs);

multispec_image = [];
for i = 1:n_imgs
    this_channel = imread(fullfile(im_dir, tifs(i).name));
    multispec_image(:,:,i) = this_channel(:,:,1);
end

rgb_image = [];
count = 1;
for i = [5,3,2]
    this_channel = imread(fullfile(im_dir, tifs(i).name));
    rgb_image(:,:,count) = this_channel(:,:,1);
    count = count + 1;
end

% mask = imread(fullfile(im_dir, tifs(end).name));
% mask = mask(:,:,1);
% instance_mask = mask == 255;

clear this_channel
% convert the images to doubles 
rgb_image = im2double(uint8(rgb_image));
multispec_image = im2double(uint8(multispec_image));
%% stretch each channel
stretched_multispec = zeros(size(multispec_image));
for i = 1:size(multispec_image,3)
    stretched_multispec(:,:,i) = imadjust(multispec_image(:,:,i));
end

%% for p28 chips, need to clip into individual chips
n_chips = 4;
square_size = 600;
multispec_chips = zeros(square_size,square_size,n_imgs,n_chips);

for i = 1:n_chips
    % display which chip to click
    to_disp = 'click the upper left of chip %d';
    str = sprintf(to_disp,i);
    disp(str)
    % get the coordinates through ginput;
    imagesc(stretched_multispec(:,:,3))
    colormap gray
    [col_coord,row_coord] = ginput(1);
    row_coord = round(row_coord);
    col_coord = round(col_coord);
    
    multispec_chips(:,:,:,i) = stretched_multispec(row_coord:(row_coord+square_size-1),col_coord:(col_coord+square_size-1),:);
end

%% pca analysis of mean color values
mean_color_mat = squeeze(mean(multispec_chips,[1,2]));
mean_color_mat = mean_color_mat(:,[1,4,2,3]);
[pcs,scores,eigens] = pca(mean_color_mat');
[rgb_pcs,rgb_scores,rgb_eigens] = pca(mean_color_mat([2,3,5],:)');

colors = get(0, 'DefaultAxesColorOrder');
figure();
for i = 1:size(scores,1)
   scatter(rgb_scores(i,1),rgb_scores(i,2),100,colors(i,:),'filled')
   hold on
   scatter(scores(i,1),scores(i,2),100,colors(i,:),'filled','d')
end

xlabel('PC1')
ylabel('PC2')
%%
chips_for_fig = [1,4,2,3];

colors = get(0, 'DefaultAxesColorOrder');

figure()
for i = 1:length(chips_for_fig)
    for j = 1:size(multispec_chips,3)
    this_chip_channel = multispec_chips(:,:,j,chips_for_fig(i));
    prctile5 = prctile(this_chip_channel,5,'all');
    prctile95 = prctile(this_chip_channel,95,'all');
    plot([prctile5,prctile95],[j-0.25+(0.1*i),j-0.25+(0.1*i)],'LineWidth',5,'Color',colors(i,:))
    hold on
    scatter(mean(this_chip_channel,'all'),j-0.25+(0.1*i),100,colors(i,:),'Filled')
    end
end
xlabel('normalized pixel intensity')
ylabel('wavelength [nm]')

figure()
for i = 1:length(chips_for_fig)
    subplot(4,2,(2*i)-1)
    imshow(multispec_chips(:,:,[5,3,2],chips_for_fig(i)))
    colormap gray
    axis image
    title(num2str(i))
    subplot(4,2,(2*i))
    imshow(multispec_chips(:,:,[7,4,1],chips_for_fig(i)))
    colormap gray
    axis image
    title(num2str(i))
end

figure()
plot(squeeze(mean(multispec_chips(:,:,:,4),[1,2]))-squeeze(mean(multispec_chips(:,:,:,1),[1,2])),[1,2,3,4,5,6,7,8])
hold on
plot(squeeze(mean(multispec_chips(:,:,:,2),[1,2]))-squeeze(mean(multispec_chips(:,:,:,3),[1,2])),[1,2,3,4,5,6,7,8])
hold on
plot(squeeze(mean(multispec_chips(:,:,:,4),[1,2]))-squeeze(mean(multispec_chips(:,:,:,3),[1,2])),[1,2,3,4,5,6,7,8])
%% make a figure or two
figure()
subplot(1,3,1)

scatter(nir_scores(:,1),nir_scores(:,2))
text(nir_scores(:,1),nir_scores(:,2),cellstr(num2str([1:n_chips]')),'VerticalAlignment',...
    'bottom','HorizontalAlignment','right')
subplot(1,3,2)
scatter(rgb_scores(:,1),rgb_scores(:,2))
text(rgb_scores(:,1),rgb_scores(:,2),cellstr(num2str([1:n_chips]')),'VerticalAlignment',...
    'bottom','HorizontalAlignment','right')
subplot(1,3,3)
scatter(multispec_scores(:,1),multispec_scores(:,2))
text(multispec_scores(:,1),multispec_scores(:,2),cellstr(num2str([1:n_chips]')),'VerticalAlignment',...
    'bottom','HorizontalAlignment','right')

% 
% figure()
% for i = 1:size(multispec_chips,4)
%    subplot(3,4,i)
%    imshow(multispec_chips(:,:,[6,7,8],i))
%    title(num2str(i))
% end
% 
% figure()
% for i = 1:size(multispec_chips,4)
%    subplot(3,4,i)
%    imshow(multispec_chips(:,:,[5,3,2],i))
%    title(num2str(i))
% end
%%
figure()
scatter(multispec_norm(:,5),multispec_norm(:,42))
text(multispec_norm(:,5),multispec_norm(:,42),cellstr(num2str([1:n_chips]')),'VerticalAlignment',...
    'bottom','HorizontalAlignment','right')
xlabel('UV contrast')
ylabel('yellow contrast')

figure()
subplot(2,3,1)
imshow(stretched_chips(:,:,[5,3,2],8))
title('chip 8')
subplot(2,3,4)
imshow(stretched_chips(:,:,[5,3,2],9))
title('chip 9')
subplot(2,3,2)
imshow(stretched_chips(:,:,[5,3,2],89))
title('chip 89')
subplot(2,3,5)
imshow(stretched_chips(:,:,[5,3,2],87))
title('chip 87')
subplot(2,3,3)
imshow(stretched_chips(:,:,[5,3,2],60))
title('chip 60')
subplot(2,3,6)
imshow(stretched_chips(:,:,[5,3,2],36))
title('chip 36')

%%
figure()
for i = 1:n_chips
    subplot(10,11,i)
    imshow(stretched_chips(:,:,[5,3,2],i))
    title(num2str(i))
end
figure()
for i = 1:n_chips
    subplot(10,11,i)
    imshow(stretched_chips(:,:,7,i))
    title(num2str(i))
end
