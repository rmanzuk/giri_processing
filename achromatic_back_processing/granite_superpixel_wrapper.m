git % Script to do superpixel oversegmentation of granite images and segment
% based upon those superpixels. Not effective, but could be used for other
% images
%
% Written 12/01/2021 R. A. Manzuk
% Last edited 12/01/2021
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
neighborhood_rad = 1;
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

%% make some random crops that will make it easier to gather superpixel training data
n_crops = 5;
[crops,crop_coords] = random_crop(rgb_im(:,:,1:3),n_crops,[2000,3000]);

%% and label superpixels from those crops
n_classes = 4;
pix_per_super = 3000;
[superpix,labeled_pix] = superpixel_trainer(crops,n_classes,pix_per_super);
%% need to get the crops from the rgb im and rotation_im
crops_rgb = struct;
crops_rgb.crop1 = rgb_im(crop_coords(1,3):crop_coords(1,4)-1, crop_coords(1,1):crop_coords(1,2)-1,:);
crops_rgb.crop2 = rgb_im(crop_coords(2,3):crop_coords(2,4)-1, crop_coords(2,1):crop_coords(2,2)-1,:);
crops_rgb.crop3 = rgb_im(crop_coords(3,3):crop_coords(3,4)-1, crop_coords(3,1):crop_coords(3,2)-1,:);
crops_rgb.crop4 = rgb_im(crop_coords(4,3):crop_coords(4,4)-1, crop_coords(4,1):crop_coords(4,2)-1,:);
crops_rgb.crop5 = rgb_im(crop_coords(5,3):crop_coords(5,4)-1, crop_coords(5,1):crop_coords(5,2)-1,:);

crops_rot = struct;
crops_rot.crop1 = all_rot_im(crop_coords(1,3):crop_coords(1,4)-1, crop_coords(1,1):crop_coords(1,2)-1,:);
crops_rot.crop2 = all_rot_im(crop_coords(2,3):crop_coords(2,4)-1, crop_coords(2,1):crop_coords(2,2)-1,:);
crops_rot.crop3 = all_rot_im(crop_coords(3,3):crop_coords(3,4)-1, crop_coords(3,1):crop_coords(3,2)-1,:);
crops_rot.crop4 = all_rot_im(crop_coords(4,3):crop_coords(4,4)-1, crop_coords(4,1):crop_coords(4,2)-1,:);
crops_rot.crop5 = all_rot_im(crop_coords(5,3):crop_coords(5,4)-1, crop_coords(5,1):crop_coords(5,2)-1,:);
%% grab the training data from those super pixels
[rgb_training_data] = get_superpix_training(crops_rgb,superpix,labeled_pix);
[rot_training_data] = get_superpix_training(crops_rot,superpix,labeled_pix);

%% train svms on those data
svm_model_rgb = fitcecoc(rgb_training_data(:,1:end-1),rgb_training_data(:,end));
svm_model_rot = fitcecoc(rot_training_data(:,1:end-1),rot_training_data(:,end));

%% train random forests on those data
n_trees = 100;
trained_forest_rgb = TreeBagger(n_trees,rgb_training_data(:,1:end-1),rgb_training_data(:,end));
trained_forest_rot = TreeBagger(n_trees,rot_training_data(:,1:end-1),rot_training_data(:,end));

%% make superpixels of the whole image
n_supers = round((size(rgb_im,1) * size(rgb_im,2))/pix_per_super);
entire_superpix = superpixels(rgb_im(:,:,1:3),n_supers);

%% use svms to classify entire image
classified_im_rgb = superpix_classify_im(trained_forest_rgb,rgb_im,entire_superpix);

classified_im_rot = superpix_classify_im(svm_model_rot,all_rot_im,entire_superpix);
