% Script to look at segmentations for the granite images in cross-polarized
% light
%
% Written 11/15/2021 R. A. Manzuk
% Last edited 11/30/2021
%% set up reading of the images
input_dir = '/Users/ryan/Dropbox (Princeton)/achromatic_project/light_table_granite_ims';
file_pattern = fullfile(input_dir, '*.tif');
tifs = dir(file_pattern);
% give the indices of each color channel
blue_inds = [1:5:26];
green_inds = blue_inds + 1;
red_inds = green_inds + 1;
nir1_inds = red_inds + 1;
nir2_inds = nir1_inds + 1;

%% and just make an rgb image with added channels for gradient magnitude and direction 
% read the channels
red_im = im2double(imread(fullfile(input_dir, tifs(red_inds(1)).name)));
green_im = im2double(imread(fullfile(input_dir, tifs(green_inds(1)).name)));
blue_im = im2double(imread(fullfile(input_dir, tifs(blue_inds(1)).name)));
% use function to make rgb image, stretched
rgb_im = concat_norm_im(red_im,green_im,blue_im);

clear red_im
clear green_im
clear blue_im

% add the gradient channels, looked like radius of 10 might be best
neighborhood_rad = 5;
[Gmag_neigh_rgb,Gdir_neigh_rgb] = neighborhood_gradients(rgb_im,neighborhood_rad);

rgb_im(:,:,end+1) = Gmag_neigh_rgb;
rgb_im(:,:,end+1) = Gdir_neigh_rgb;


%% assemble a false color image with all of the rotations that includes gradient channels

% because of small channel offsets, we need to keep track of size for
% concatenation
goal_size = size(rgb_im,1);

% loop through each rotation and make an rgb, grab the gradients,
% concatenate

for i =  1:numel(green_inds)-1
    tic
    % load channels
    green_im = im2double(imread(fullfile(input_dir, tifs(green_inds(i)).name)));
    blue_im = im2double(imread(fullfile(input_dir, tifs(blue_inds(i)).name)));
    red_im = im2double(imread(fullfile(input_dir, tifs(red_inds(i)).name)));
    
    % concat and normalize the rgb for this rotation
    red_green_blue = concat_norm_im(red_im([end-goal_size+1]:end,:),green_im([end-goal_size+1]:end,:),blue_im([end-goal_size+1]:end,:));

    % get gradient stuff
    [Gmag_neigh_rgb,Gdir_neigh_rgb] = neighborhood_gradients(red_green_blue, neighborhood_rad);

    % make 5 channel
    red_green_blue(:,:,end+1) = Gmag_neigh_rgb;
    red_green_blue(:,:,end+1) = Gdir_neigh_rgb;
    
    % concatenate into the big multichannel conglomerate
    if i == 1
        all_rot_im = red_green_blue;
    else
        all_rot_im = cat(3,all_rot_im,red_green_blue);
    end
    
    clear red_green_blue
    clear red_im
    clear blue_im
    clear green_im
    toc
end

%% ginput the training boxes
% clicking training data in the order quartz, orthoclase, plagioclase,
% mafics.
[~, rectangle_coords] = extract_image_rectangles(rgb_im(:,:,1:3),4);
%% grab those training data for the rgb and rotation image
[rgb_training_data, ~] = extract_image_rectangles(rgb_im,4, rectangle_coords);
[rot_training_data, ~] = extract_image_rectangles(all_rot_im,4, rectangle_coords);

%% (optional) add to those training rectangles
[pixel_vals2,rectangle_coords2] = add_image_rectangles(rgb_im,rgb_training_data,rectangle_coords);
%% train a multisvm on the rgb data and classify
rgb_column = reshape(rgb_im,size(rgb_im,1)*size(rgb_im,2),size(rgb_im,3));
svm_model_rgb = fitcecoc(rgb_training_data(:,1:end-1),rgb_training_data(:,end));

predicted_rgb_image = predict(svm_model_rgb,rgb_column);
class_image_rgb = reshape(predicted_rgb_image,size(rgb_im,1),size(rgb_im,2));

%% train a multisvm on all the rotated channels and classify
rot_im_column = reshape(all_rot_im,size(all_rot_im,1)*size(all_rot_im,2),size(all_rot_im,3));
svm_model_rot = fitcecoc(rot_training_data(:,1:end-1),rot_training_data(:,end));

predicted_rot_image = predict(svm_model_rot,rot_im_column);
class_image_rot = reshape(predicted_rot_image,size(all_rot_im,1),size(all_rot_im,2));

%% perform morphological operations on classified images
rad = 2;
n_pix = 1500;
class_image_rot2 = multiclass_open(class_image_rot,rad,n_pix);
class_image_rgb2 = multiclass_open(class_image_rgb,n_pix);
%% make a nice little figure
figure;
subplot(2,1,1)
imshow(rgb_im(:,:,[1:3]))
%plot_class_rectangles(rgb_im,rectangle_coords)
axis image
subplot(2,2,3)
imagesc(class_image_rgb)
title('Single RGB')
axis image
subplot(2,2,4)
imagesc(class_image_rot)
title('All rotations RGB')
axis image


%% traced image for accuracy 
% testing image coordinates are 4692:6180, 9594:12906
% indeces: quartz=4, orth=5, plag=3,lith=2
%load it
ground_truth = imread('/Users/ryan/Dropbox (Princeton)/achromatic_project/light_table_granite_ims/photoshop_tracing/testing_im_smaller.tif');

% switch the indeces to match what we've been doing
ground_truth(ground_truth == 4) = 1;
ground_truth(ground_truth == 2) = 4;
ground_truth(ground_truth == 5) = 2;

% set some stuff to NaN because it wasn't traced
ground_truth(ground_truth == 5) = NaN;
ground_truth(ground_truth == 6) = NaN;
ground_truth(ground_truth == 0) = NaN;
% how many pixels were actually traced?
n_traced = sum(and(ground_truth > 0, ground_truth < 5),'all');

% and compare segmentations to test tracing for accuracy
rgb_test_seg = class_image_rgb(4692:6180,9594:12906);
rot_test_seg = class_image_rot2(4692:6180,9594:12906);

rgb_acc = (sum(rgb_test_seg == ground_truth,'all')/n_traced) * 100;
rot_acc = (sum(rot_test_seg == ground_truth,'all')/n_traced) * 100;

%% what is the impact of field of view and point counting?
% fov of a microscope camera is 2048x2048
fov_size = [2048,2048];
n_samples = 100;
n_points = 100;
[class_seg_fractions_2_5, class_stoch_fracs_2_5, class_fixed_fracs_2_5] = fov_sample(class_image_rot, fov_size, n_samples,'PointCount','both',n_points);
[class_seg_fractions_10, class_stoch_fracs_10, class_fixed_fracs_10] = fov_sample(class_image_rot, fov_size/4, n_samples,'PointCount','both',n_points);
n_classes = 4;

% what are the real fractions in the segmentation?
real_fracs = zeros(1,n_classes);
for j = 1:n_classes
    pix_present = sum(class_image_rot == j, 'all');
    real_fracs(1,j) = pix_present / numel(class_image_rot);
end

%% make a figure on impact of field of view and point counting
% jitter plots of the different methods of assessing modality
figure();
labels = [1;2*ones(n_samples,1);3*ones(n_samples,1);4*ones(n_samples,1);5*ones(n_samples,1);...
    6;7*ones(n_samples,1);8*ones(n_samples,1);9*ones(n_samples,1);10*ones(n_samples,1);...
    11;12*ones(n_samples,1);13*ones(n_samples,1);14*ones(n_samples,1);15*ones(n_samples,1);...
    16;17*ones(n_samples,1);18*ones(n_samples,1);19*ones(n_samples,1);20*ones(n_samples,1)];
to_jitter = [real_fracs(1); class_seg_fractions_2_5(:,1); class_fixed_fracs_2_5(:,1);...
    class_seg_fractions_10(:,1); class_fixed_fracs_10(:,1);real_fracs(2); class_seg_fractions_2_5(:,2); class_fixed_fracs_2_5(:,2);...
    class_seg_fractions_10(:,2); class_fixed_fracs_10(:,2);real_fracs(3); class_seg_fractions_2_5(:,3); class_fixed_fracs_2_5(:,3);...
    class_seg_fractions_10(:,3); class_fixed_fracs_10(:,3);real_fracs(4); class_seg_fractions_2_5(:,4); class_fixed_fracs_2_5(:,4);...
    class_seg_fractions_10(:,4); class_fixed_fracs_10(:,4);];
jitter = 0.2;
colors = [repmat([255/255,255/255,255/255],5,1);repmat([10/255,107/255,108/255],5,1);...
    repmat([145/255,167/255,126/255],5,1);repmat([32/255,25/255,55/255],5,1)];
groups = repmat([1,2,3,2,3],1,4);
jitterplot(to_jitter,labels,jitter,colors,groups);