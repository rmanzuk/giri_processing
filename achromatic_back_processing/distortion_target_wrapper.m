% script to run through functions and get data from the distortion target,
% specifically checking if there is a distortion field when switching
% between wavelengths.
%
% R. A. Manzuk 01/07/2021
%% load everything

focus_im = imread('/Users/ryan/Dropbox (Princeton)/achromatic_project/petro_scale_photos/Image0000.tif');
%focus_im = rgb2gray(focus_im);
focus_im = im2double(focus_im);

blur_im = imread('/Users/rmanzuk/Desktop/achromatic_project/test_images/625p_red_focus_target6019_1.tif');
%blur_im = rgb2gray(blur_im);
blur_im = im2double(blur_im);

%% dots are blobs, lets blob detect
% for nonmaximum suppression
local_thresh = 3;
% size of the blobs
rad = 50;
% how strong of a response defines a blob
global_thresh = 0.05;

focus_blobs = detect_blobs(focus_im,rad,local_thresh,global_thresh);
blur_blobs = detect_blobs(blur_im,rad,local_thresh,global_thresh);
% and bring blobs from image space to xy space
focus_blobs_xy = [focus_blobs(:,2),-focus_blobs(:,1)];
blur_blobs_xy = [blur_blobs(:,2),-blur_blobs(:,1)];
%% find the nearest neighbors
all_distances = pdist2(focus_blobs_xy,blur_blobs_xy);
[~,closest_ind] = min(all_distances);
clear all_distances
matched_focus_blobs = focus_blobs_xy(closest_ind',:);

%% make into a vector field
distortion_vectors = [blur_blobs_xy(:,1)-matched_focus_blobs(:,1),blur_blobs_xy(:,2)-matched_focus_blobs(:,2)];
vector_magnitudes = sqrt((distortion_vectors(:,1).^2) + (distortion_vectors(:,2).^2));
okay_ind = vector_magnitudes <50;
figure(1);
subplot(2,1,1)
quiver(matched_focus_blobs(okay_ind,1),matched_focus_blobs(okay_ind,2),distortion_vectors(okay_ind,1),distortion_vectors(okay_ind,2));
axis tight
subplot(2,1,2)
scatter(matched_focus_blobs(okay_ind,1),matched_focus_blobs(okay_ind,2),20,vector_magnitudes(okay_ind),'Filled');
colorbar
axis tight
%% what kernel makes our blurred image?
% images are too big to process at once...tile them
[blur_tiles,blur_centers] = tile_image(blur_im,1000);
[focus_tiles,focus_centers] = tile_image(focus_im,1000);

% how big of a kernel to estimate (arbitrary for now)
kernel_size = [3,3];

% empty matrix to accept the kernel for each tile
final_kernels = zeros(kernel_size(1),kernel_size(2),size(blur_tiles,3));
for i = 1:size(blur_tiles,3)
    tic
    % which tiles are we working with
    this_blur_tile = blur_tiles(:,:,i);
    this_focus_tile = focus_tiles(:,:,i);
    
    % convolution removes some info from edges, so blurred tiles need outer
    % pixels removed
    blurred_valid = this_blur_tile(ceil(kernel_size(1)/2):end-floor(kernel_size(1)/2), ...
          ceil(kernel_size(2)/2):end-floor(kernel_size(2)/2));
    
    % make a column for each chunk of the sharp image that is input to the
    % blurring convolution
    sharp_blocks = im2col(this_focus_tile, kernel_size, 'sliding')';
    
    % We've now set up a least squares problem...solve it
    xcor_kernel = sharp_blocks \ blurred_valid(:);
    
    % make the kernel from a column into a 2d matrix
    final_kernels(:,:,i)= rot90(reshape(xcor_kernel, kernel_size), 2);
    toc
end