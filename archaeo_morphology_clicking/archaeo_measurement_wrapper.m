% script to run through different blocks/functions to analyze different archaeos

% R. A. Manzuk 09/18/2020
%% Load necessary stuff for Stuart's Mill Sample
load('sm_all_inners.mat')
load('sm_all_outers.mat')
block_top_sd = [35,10];
strike_im_heading = 307;
imageDir = dir(fullfile('/Users/rmanzuk/Desktop/nevada_jan_2019/sm_117_71_downsampled_25/', '*.tif'));
bedding_sd = [238, 34];
scale_ratio = 6.874/2;

inners = sm_inners;
outers = sm_outers;

%% Load necessary stuff for RLG136a
load('rlg136a_all_inners.mat')
load('rlg136a_all_outers.mat')
scale_ratio = 6.874;
block_top_sd = [194,63];
strike_im_heading = 333;
bedding_sd = [298, 23];
imageDir = dir(fullfile('/Users/rmanzuk/Desktop/ketza_2019/grinder_stacks/rlg_136a/downsampled_stack', '*.tif'));

inners = rlg136a_inners;
outers = rlg136a_outers;
%% Before we do anything, let's densify slices and get the branching points
[inner_dense_slices,outer_dense_slices] = densify_slices(inners,outers,3);
[branched_flags,branching_angles,branch_points_3d] = process_branched_network(sm_inners,sm_outers,scale_ratio,10);

%% Now we can densify in 3d, and take the center lines
[inner_3d_dense,outer_3d_dense] = densify_3d(inner_dense_slices,outer_dense_slices,scale_ratio,3,1);
[inner_center_points,outer_center_points] = easy_center_lines(inner_3d_dense,outer_3d_dense,10,50,1);

%% and rotate everybody
[inner_centers_rotated, outer_centers_rotated] = rotate_clicked_data3d(inner_center_points, outer_center_points, block_top_sd, strike_im_heading, bedding_sd,0);
[inner_3d_rotated, outer_3d_rotated] = rotate_clicked_data3d(inner_3d_dense, outer_3d_dense, block_top_sd, strike_im_heading, bedding_sd,1);
[branch_points_rotated] = rotate_branch_points(branch_points_3d, block_top_sd, strike_im_heading, bedding_sd);

%% and get some data
sampling_freq = 10;
thickness_samp = 1;
[inner_center_stats,outer_center_stats] = center_line_analysis(inner_3d_rotated,outer_3d_rotated,outer_centers_rotated,outer_centers_rotated,sampling_freq,thickness_samp);