% script to run through different blocks/functions to analyze different archaeos

% R. A. Manzuk 09/18/2020
% Last edited: 03/01/2021
%% Load necessary stuff for Stuart's Mill Sample
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/sm_all_inners.mat')
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/sm_all_outers.mat')
sm_block_top_sd = [35,10];
sm_strike_im_heading = 307;
sm_bedding_sd = [238, 34];
sm_scale_ratio = 6.874/2;
sm_um_pixel = 145.4;

%% Load necessary stuff for cc297
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/cc297_all_inners.mat')
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/cc297_all_outers.mat')
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/cc297_crack_left.mat')
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/cc297_crack_right.mat')

cc297_scale_ratio = 10;
cc297_block_top_sd = [180,90];
cc297_strike_im_heading = 90;
cc297_bedding_sd = [170, 10];
cc297_um_pixel = 40.4;
%% Load necessary stuff for Labrador Sample
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/labrador_all_inners.mat');
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/labrador_all_outers.mat');

labrador_block_top_sd = [127,94];
labrador_strike_im_heading = 68;
labrador_bedding_sd = [187, 15];
labrador_scale_ratio = 6.874;
labrador_um_pixel = 72.7;

%% Run through all the parts for stewart's mill
% densify the slices
tic
sm_inner_dense_slices = densify_slices(sm_inners,5);
sm_outer_dense_slices = densify_slices(sm_outers,5);
disp('done densifying slices')
toc

% identify branch points
tic
[sm_branched_flags,sm_branch_points_3d] = process_branched_network(sm_outers,sm_scale_ratio);
disp('done identifying branch points')
toc

% densify and make 3d
points_here_thresh = 20;
iterate_stop = 4;
sample_freq = 5;
thickness_sampling = 20;
n_iter = 1;
densification = 2;

tic
sm_inner_3d_dense = densify_3d(sm_inner_dense_slices,sm_scale_ratio,densification);
sm_outer_3d_dense = densify_3d(sm_outer_dense_slices,sm_scale_ratio,densification);
disp('done making 3d and densifying')
toc

% take the center lines
tic
sm_center_points = iterate_center_lines(sm_outer_3d_dense,points_here_thresh,iterate_stop,n_iter,sample_freq,thickness_sampling,15);
sm_center_points = cull_centers(sm_center_points,1);
disp('done taking the center lines')
toc

% rotate everybody
tic
sm_centers_rotated = rotate_clicked_data3d(sm_center_points, sm_block_top_sd, sm_strike_im_heading, sm_bedding_sd);
sm_inner_3d_rotated = rotate_clicked_data3d(sm_inner_3d_dense, sm_block_top_sd, sm_strike_im_heading, sm_bedding_sd);
sm_outer_3d_rotated = rotate_clicked_data3d(sm_outer_3d_dense, sm_block_top_sd, sm_strike_im_heading, sm_bedding_sd);
sm_branch_points_rotated = rotate_branch_points(sm_branch_points_3d, sm_block_top_sd, sm_strike_im_heading, sm_bedding_sd);
disp('done rotating data')
toc

% and get some data
tic
sampling_freq = 5;
thickness_samp = 4;
diff_thresh = 5;

sm_inner_center_stats = center_line_analysis(sm_inner_3d_rotated,sm_centers_rotated,sample_freq,thickness_sampling);
sm_outer_center_stats = center_line_analysis(sm_outer_3d_rotated,sm_centers_rotated,sample_freq,thickness_sampling);

[sm_deriv_means,sm_deriv_variances,sm_thicks_encountered,sm_nn_dists] = centers_plane_pass(sm_outer_center_stats,sm_outer_3d_rotated, diff_thresh);
[sm_mean_declinations,sm_mean_slope_runs] = cart2pol(sm_deriv_means(:,1),sm_deriv_means(:,2));
sm_mean_inclinations = atand(sm_deriv_means(:,3)./sm_mean_slope_runs);

[sm_br_angles, sm_br_lengths] = spline_branch_angles2(sm_branch_points_rotated,sm_branched_flags,sm_outer_center_stats, 10);

[sm_surface_area,sm_volume] = sa_and_vol(sm_outer_dense_slices,sm_scale_ratio,sm_um_pixel/1e4,sm_branched_flags);
[sm_footprint_points, sm_footprint_area] = get_footprint(sm_outer_3d_rotated, 3, sm_um_pixel/1e4);
[sm_convhull_points, sm_enclosing_volume] = get_enclosing_volume(sm_outer_3d_rotated, sm_um_pixel/1e4);
[sm_surf_areas] = std_sa_encvol(sm_outer_dense_slices,sm_scale_ratio,sm_um_pixel,10,30000,sm_branched_flags);
disp('done gathering data')
toc

%% Run through all the parts for ketza cc297
% densify the slices
tic
cc297_inner_dense_slices = densify_slices(cc297_inners,5);
cc297_outer_dense_slices = densify_slices(cc297_outers,5);
disp('done densifying slices')
toc

% identify branch points
tic
[cc297_branched_flags,cc297_branch_points_3d] = process_branched_network(cc297_outers,cc297_scale_ratio);
disp('done identifying branch points')
toc

% densify and make 3d
points_here_thresh = 20;
iterate_stop = 4;
sample_freq = 5;
thickness_sampling = 30;
n_iter = 1;
densification = 2;

tic
cc297_inner_3d_dense = densify_3d(cc297_inner_dense_slices,cc297_scale_ratio,densification);
cc297_outer_3d_dense = densify_3d(cc297_outer_dense_slices,cc297_scale_ratio,densification);
crack_left_3d = densify_3d(cc297_crack_left,cc297_scale_ratio,densification);
crack_right_3d = densify_3d(cc297_crack_right,cc297_scale_ratio,densification);
disp('done making 3d and densifying')
toc

tic
[cc297_collapsed_outers, cc297_collapsed_inners, cc297_collapsed_branch_points] = collapse_calcite(crack_left_3d, crack_right_3d, cc297_outer_3d_dense, cc297_inner_3d_dense, cc297_branch_points_3d);
disp('done collapsing calcite')
toc 

% take the center lines
tic
sampling_resolution = 30;
cc297_center_points = iterate_center_lines(cc297_collapsed_outers,points_here_thresh,iterate_stop,n_iter,sample_freq,thickness_sampling,20);
cc297_center_points = cull_centers(cc297_center_points,1);
disp('done taking the center lines')
toc

% rotate everybody
tic
cc297_centers_rotated = rotate_clicked_data3d(cc297_center_points, cc297_block_top_sd, cc297_strike_im_heading, cc297_bedding_sd);
cc297_inner_3d_rotated = rotate_clicked_data3d(cc297_collapsed_inners, cc297_block_top_sd, cc297_strike_im_heading, cc297_bedding_sd);
cc297_outer_3d_rotated = rotate_clicked_data3d(cc297_collapsed_outers, cc297_block_top_sd, cc297_strike_im_heading, cc297_bedding_sd);
cc297_branch_points_rotated = rotate_branch_points(cc297_branch_points_3d, cc297_block_top_sd, cc297_strike_im_heading, cc297_bedding_sd);
disp('done rotating data')
toc

% and get some data
tic
sampling_freq = 5;
thickness_samp = 4;
diff_thresh = 5;

cc297_inner_center_stats = center_line_analysis(cc297_inner_3d_rotated,cc297_centers_rotated,sample_freq,thickness_sampling);
cc297_outer_center_stats = center_line_analysis(cc297_outer_3d_rotated,cc297_centers_rotated,sample_freq,thickness_sampling);

[cc297_deriv_means,cc297_deriv_variances,cc297_thicks_encountered,cc297_nn_dists] = centers_plane_pass(cc297_outer_center_stats,cc297_outer_3d_rotated, diff_thresh);
[cc297_mean_declinations,cc297_mean_slope_runs] = cart2pol(cc297_deriv_means(:,1),cc297_deriv_means(:,2));
cc297_mean_inclinations = atand(cc297_deriv_means(:,3)./cc297_mean_slope_runs);

[cc297_br_angles, cc297_br_lengths]  = spline_branch_angles2(cc297_branch_points_rotated,cc297_branched_flags,cc297_outer_center_stats, 10);

[cc297_surface_area,cc297_volume] = sa_and_vol(cc297_outer_dense_slices,cc297_scale_ratio,cc297_um_pixel/1e4,cc297_branched_flags);
[cc297_footprint_points, cc297_footprint_area] = get_footprint(cc297_outer_3d_rotated, 3, cc297_um_pixel/1e4);
[cc297_convhull_points, cc297_enclosing_volume] = get_enclosing_volume(cc297_outer_3d_rotated, cc297_um_pixel/1e4);
[cc297_surf_areas] = std_sa_encvol(cc297_outer_dense_slices,cc297_scale_ratio,cc297_um_pixel,10,30000,cc297_branched_flags);
disp('done gathering data')
toc

%% Run through all the parts for labrador
% densify the slices
tic
labrador_inner_dense_slices = densify_slices(labrador_inners,5);
labrador_outer_dense_slices = densify_slices(labrador_outers,5);
disp('done densifying slices')
toc

% identify branch points
tic
[labrador_branched_flags,labrador_branch_points_3d] = process_branched_network(labrador_outers,labrador_scale_ratio);
disp('done identifying branch points')
toc

% densify and make 3d
points_here_thresh = 20;
iterate_stop = 4;
sample_freq = 5;
thickness_sampling = 40;
n_iter = 1;
densification = 2;

tic
labrador_inner_3d_dense = densify_3d(labrador_inner_dense_slices,labrador_scale_ratio,densification);
labrador_outer_3d_dense = densify_3d(labrador_outer_dense_slices,labrador_scale_ratio,densification);
disp('done making 3d and densifying')
toc

% take the center lines
tic

labrador_center_points = iterate_center_lines(labrador_outer_3d_dense,points_here_thresh,iterate_stop,n_iter,sample_freq,thickness_sampling,50);
labrador_center_points = cull_centers(labrador_center_points,0.3);
disp('done taking the center lines')
toc

% rotate everybody
tic
labrador_centers_rotated = rotate_clicked_data3d(labrador_center_points, labrador_block_top_sd, labrador_strike_im_heading, labrador_bedding_sd);
labrador_inner_3d_rotated = rotate_clicked_data3d(labrador_inner_3d_dense, labrador_block_top_sd, labrador_strike_im_heading, labrador_bedding_sd);
labrador_outer_3d_rotated = rotate_clicked_data3d(labrador_outer_3d_dense, labrador_block_top_sd, labrador_strike_im_heading, labrador_bedding_sd);
labrador_branch_points_rotated = rotate_branch_points(labrador_branch_points_3d, labrador_block_top_sd, labrador_strike_im_heading, labrador_bedding_sd);
disp('done rotating data')
toc

% and get some data
tic
sampling_freq = 5;
thickness_samp = 4;
diff_thresh = 5;

labrador_inner_center_stats = center_line_analysis(labrador_inner_3d_rotated,labrador_centers_rotated,sample_freq,thickness_sampling);
labrador_outer_center_stats = center_line_analysis(labrador_outer_3d_rotated,labrador_centers_rotated,sample_freq,thickness_sampling);

[labrador_deriv_means,labrador_deriv_variances,labrador_thicks_encountered,labrador_nn_dists] = centers_plane_pass(labrador_outer_center_stats,labrador_outer_3d_rotated, diff_thresh);
[labrador_mean_declinations,labrador_mean_slope_runs] = cart2pol(labrador_deriv_means(:,1),labrador_deriv_means(:,2));
labrador_mean_inclinations = atand(labrador_deriv_means(:,3)./labrador_mean_slope_runs);

[labrador_br_angles, labrador_br_lengths] = spline_branch_angles2(labrador_branch_points_rotated,labrador_branched_flags,labrador_outer_center_stats, 10);

[labrador_surface_area,labrador_volume] = sa_and_vol(labrador_outer_dense_slices,labrador_scale_ratio,labrador_um_pixel/1e4,labrador_branched_flags);
[labrador_footprint_points, labrador_footprint_area] = get_footprint(labrador_outer_3d_rotated, 3, labrador_um_pixel/1e4);
[labrador_convhull_points, labrador_enclosing_volume] = get_enclosing_volume(labrador_outer_3d_rotated, labrador_um_pixel/1e4);
[labrador_surf_areas] = std_sa_encvol(labrador_outer_dense_slices,labrador_scale_ratio,labrador_um_pixel,10,30000,labrador_branched_flags);
disp('done gathering data')
toc

