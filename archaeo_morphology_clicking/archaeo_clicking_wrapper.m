% script with all of the necessary components to load and click any of the
% archaeo samples

% R. A. Manzuk 09/18/2020
% Last edited: 03/01/2021
%% Load necessary stuff for Stuart's Mill Sample
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/sm_all_inners.mat')
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/sm_all_outers.mat')
block_top_sd = [35,10];
strike_im_heading = 307;
input_folder = '/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/grinder_stacks/sm_117_71_downsampled_25';
bedding_sd = [238, 34];
scale_ratio = 6.874/2;
um_pixel = 145.4;


%% Load necessary stuff for cc297
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/cc297_all_inners.mat')
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/cc297_all_outers.mat')
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/cc297_crack_left.mat')
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/cc297_crack_right.mat')

scale_ratio = 10;
block_top_sd = [0,90];
strike_im_heading = 270;
bedding_sd = [193, 10];
input_folder = '/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/grinder_stacks/cc297/transverse_resample_every10';
um_pixel = 40.4;

%% Load necessary stuff for Labrador Sample
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/labrador_all_inners.mat');
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/labrador_all_outers.mat');

block_top_sd = [127,94];
strike_im_heading = 68;
input_folder = '/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/grinder_stacks/labrador_r02/8bit_quarter_scale_every25';
bedding_sd = [187, 15];
scale_ratio = 6.874;
um_pixel = 72.7;

%% Load necessary stuff for RLG136a (did not use this sample)
% load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/rlg136a_all_inners.mat')
% load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/rlg136a_all_outers.mat')
% scale_ratio = 6.874;
% block_top_sd = [194,63];
% strike_im_heading = 333;
% bedding_sd = [298, 23];
% input_folder = '/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/grinder_stacks/rlg_136a';
% um_pixel = 72.7;