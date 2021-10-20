% This script runs through the best_blur.m function with the proper set of
% images four our iq4 experiments to get the best sigma for gaussian
% simulated blur between each pair of wavelengths.

% Ryan A. Manzuk. Created 01/28/2021
% Last edited 02/11/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% set up image reading
input_dir = '/Users/rmanzuk/Desktop/achromatic_project/mansouri_image_set';
file_pattern = fullfile(input_dir, '*.tif');
tifs = dir(file_pattern);

n_imgs = numel(tifs);

% in this directory, images have the naming scheme
% sharpWavelength_imageWavelength. for example 530_470 would be an image
% taken at 470nm at the table height that is sharpest for 530 nm
name_sharp_lambda = zeros(1,n_imgs);
name_this_lambda = zeros(1,n_imgs);
for i = 1:n_imgs
    name_sharp_lambda(i) = str2num(tifs(i).name(1:3));
    name_this_lambda(i) = str2num(tifs(i).name(5:7));
end
%%
sigma = [0.1:0.1:10];
wavelengths = [470,530,590,625,730,850,940];
mse_cross_mat = zeros(sqrt(n_imgs),sqrt(n_imgs),length(sigma));

% conv net for denoising
net = denoisingNetwork('DnCNN');

% sharp images are i, blurry images are j
for i = 1:size(mse_cross_mat,1)
    for j = 1:size(mse_cross_mat,2)
        
        % find index in the image file list for each image.
        sharp_image_ind = (name_sharp_lambda == wavelengths(i)) & (name_this_lambda == wavelengths(i));
        blur_image_ind = (name_sharp_lambda == wavelengths(i)) & (name_this_lambda == wavelengths(j));
        
        % read the proper images
        sharp_im = imread(fullfile(input_dir, tifs(sharp_image_ind).name));
        blur_im = imread(fullfile(input_dir, tifs(blur_image_ind).name));
        
        % tile the images to make them smaller
        [sharp_tiles,~] = tile_image(sharp_im,1000);
        [blur_tiles,~] = tile_image(blur_im,1000);
        
        % go through each tile of these two images and assess the MSE for
        % all elements of sigma
        these_MSEs = zeros(size(sharp_tiles,3),length(sigma));
        for k = [1:10:size(sharp_tiles,3)]
            tic
            these_MSEs(k,:) = best_blur(sharp_tiles(:,:,k),blur_tiles(:,:,k),sigma,net);
            toc
        end
        % take the mean over all tiles
        mean_MSE = mean(these_MSEs);
        mse_cross_mat(i,j,:) = reshape(mean_MSE,1,1,length(sigma));
    end
end
%%

% plotting all of the MSE vs. sigma curves.
MSE_min_sigma = zeros(size(mse_cross_mat,1),size(mse_cross_mat,2));
for i = 1:size(mse_cross_mat,1)
    for j = 1:size(mse_cross_mat,2)
        this_MSE_signal = squeeze(mse_cross_mat(i,j,:));
        [~,min_ind] = min(this_MSE_signal);
        MSE_min_sigma(i,j) = sigma(min_ind);
        plot(sigma,this_MSE_signal)
        hold on
    end
end
