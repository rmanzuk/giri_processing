% scrip to hopefully automate the tracing of edges in coral ct scans
% (analagous to manual clicking in grinder stacks) to then put traces
% through same measurement pipeline.

% R. A. Manzuk 12/15/2020
%% set up the image directory

input_folder = '/Users/ryan/Desktop/branch_angle_project/coral_ct_data/kaandorp_scans/8865_8866_raw_tifs_resampled';

file_pattern = fullfile(input_folder, '*.tif');
tifs = dir(file_pattern);
base_names = natsortfiles({tifs.name});
fulln=numel(base_names);
%% loop though all images and extract edges

% tunable parameters
blur_kernel_size = 2;
island_size = 50;
hole_size = 500;
distance_thresh = 10;
binary_thresh = 0.2;

% set up final_outers to receive edges
final_outers = {};
% set up a counter
counter = 0;
for i = 150:415%numel(base_names)
    %read the image
    this_im = imread(fullfile(input_folder, base_names{i}));
    
    % make the image a double and stretch it a bit
    im_double = im2double(this_im)./max(im2double(this_im),[],'all');
    % use function to get this slice's edges
    edge_coords = get_slice_edges2(im_double,blur_kernel_size,island_size,hole_size,binary_thresh);
    % now, the fun, deciding which edges belong to the same branch
    
    % if there were no edges in the image, we don't need to do
    % anything...so ask if it's not empty
    if ~isempty(edge_coords)
        % we have edges here, so update the counter
        counter = counter + 1;
        if counter == 1
            % if these are our first edges, we just take them, no questions
            % asked
            final_outers{i} = edge_coords;
        else
            % otherwise we have to match the edges, use the function
            final_outers{i} = match_edges(edge_coords,final_outers{i-1},distance_thresh);
            
        end
    else
        %do nothing if there are no edges
        final_outers{i} = {};
    end
end
%% Cleaning up the automated tracing output
cleaned_outers = remove_spurious_edges(final_outers);

%% load pre-processed, cleaned data

load('/Users/ryan/Desktop/branch_angle_project/coral_ct_data/millepora_cleaned.mat');
load('/Users/ryan/Desktop/branch_angle_project/coral_ct_data/loripes_cleaned.mat');
load('/Users/ryan/Desktop/branch_angle_project/coral_ct_data/cytherea_cleaned.mat');
load('/Users/ryan/Desktop/branch_angle_project/coral_ct_data/caroliana_cleaned.mat');

load('/Users/ryan/Desktop/branch_angle_project/coral_ct_data/madracis6m_cleaned.mat');
load('/Users/ryan/Desktop/branch_angle_project/coral_ct_data/madracis15m_cleaned.mat');
load('/Users/ryan/Desktop/branch_angle_project/coral_ct_data/madracis20m_cleaned.mat');

%%
% rearrange the cell array to be 1xn_branches
millepora_reshaped = reshape_coral_cell(millepora_cleaned);
loripes_reshaped = reshape_coral_cell(loripes_cleaned);
cytherea_reshaped = reshape_coral_cell(cytherea_cleaned);
caroliana_reshaped = reshape_coral_cell(caroliana_cleaned);

madracis6m_reshaped = reshape_coral_cell(madracis6m_cleaned);
madracis15m_reshaped = reshape_coral_cell(madracis15m_cleaned);
madracis20m_reshaped = reshape_coral_cell(madracis20m_cleaned);

clear millepora_cleaned
clear loripes_cleaned
clear cytherea_cleaned
clear caroliana_cleaned
clear madracis6m_cleaned
clear madracis15m_cleaned
clear madracis20m_cleaned


% we also need the edges sorted as if going around the circle, not with
% sequential indices
millepora_resorted = sort_outline_points(millepora_reshaped);
loripes_resorted = sort_outline_points(loripes_reshaped);
cytherea_resorted = sort_outline_points(cytherea_reshaped);
caroliana_resorted = sort_outline_points(caroliana_reshaped);

madracis6m_resorted = sort_outline_points(madracis6m_reshaped);
madracis15m_resorted = sort_outline_points(madracis15m_reshaped);
madracis20m_resorted = sort_outline_points(madracis20m_reshaped);

clear millepora_reshaped
clear loripes_reshaped
clear cytherea_reshaped
clear caroliana_reshaped
clear madracis6m_reshaped
clear madracis15m_reshaped
clear madracis20m_reshaped


% and then let's downsample the slices
densify_factor = 0.3;
 
millepora_downsampled = densify_slices(millepora_resorted,0.2);
loripes_downsampled = densify_slices(loripes_resorted,densify_factor);
cytherea_downsampled = densify_slices(cytherea_resorted,densify_factor);
caroliana_downsampled = densify_slices(caroliana_resorted,densify_factor);

madracis6m_downsampled = densify_slices(madracis6m_resorted,0.8);
madracis15m_downsampled = densify_slices(madracis15m_resorted,0.8);
madracis20m_downsampled = densify_slices(madracis20m_resorted,1);
clear millepora_resorted
clear loripes_resorted
clear cytherea_resorted
clear caroliana_resorted
clear madracis6m_resorted
clear madracis15m_resorted
clear madracis20m_resorted
%% The scales of each specimen in um
% caroliniana 
caroliana_scale = 81.1170;
% cytherea
cytherea_scale = 107.8580;
% loripes 
loripes_scale = 96.5800;
% millepora 
millepora_scale = 49.9990;

% madracis 6m
madracis6m_scale = 250;
% madracis 15m
madracis15m_scale = 250;
% madracis 20m
madracis20m_scale = 250;
%% and get the branching points
distance_threshold = 20;
combining_threshold = 50;
scale_ratio = 1;
[millepora_branched_flags,millepora_branching_points_3d,millepora_combined] = id_branch_points(millepora_downsampled,scale_ratio,distance_threshold,combining_threshold);
[loripes_branched_flags,loripes_branching_points_3d,loripes_combined] = id_branch_points(loripes_downsampled,scale_ratio,distance_threshold,combining_threshold);
[cytherea_branched_flags,cytherea_branching_points_3d,cytherea_combined] = id_branch_points(cytherea_downsampled,scale_ratio,distance_threshold,combining_threshold);
[caroliana_branched_flags,caroliana_branching_points_3d,caroliana_combined] = id_branch_points(caroliana_downsampled,scale_ratio,distance_threshold,combining_threshold);

[madracis6m_branched_flags,madracis6m_branching_points_3d,madracis6m_combined] = id_branch_points(madracis6m_downsampled,scale_ratio,distance_threshold,combining_threshold);
[madracis15m_branched_flags,madracis15m_branching_points_3d,madracis15m_combined] = id_branch_points(madracis15m_downsampled,scale_ratio,distance_threshold,combining_threshold);
[madracis20m_branched_flags,madracis20m_branching_points_3d,madracis20m_combined] = id_branch_points(madracis20m_downsampled,scale_ratio,distance_threshold,combining_threshold);
% clear millepora_downsampled
% clear loripes_downsampled
% clear cytherea_downsampled
% clear caroliana_downsampled
%% take slice data and make 3d
millepora_3d = make_clicking_3d(millepora_combined,scale_ratio);
loripes_3d = make_clicking_3d(loripes_combined,scale_ratio);
cytherea_3d = make_clicking_3d(cytherea_combined,scale_ratio);
caroliana_3d = make_clicking_3d(caroliana_combined,scale_ratio);

madracis6m_3d = make_clicking_3d(madracis6m_combined,scale_ratio);
madracis15m_3d = make_clicking_3d(madracis15m_combined,scale_ratio);
madracis20m_3d = make_clicking_3d(madracis20m_combined,scale_ratio);
% clear millepora_combined
% clear loripes_combined
% clear cytherea_combined
% clear caroliana_combined
%% get the center lines
sampling_resolution = 8;
sampling_freq = 3;
points_here_thresh = 30;

millepora_center_points = easy_center_lines(millepora_3d,20,sampling_freq,points_here_thresh);
loripes_center_points = easy_center_lines(loripes_3d,20,sampling_freq,points_here_thresh);
cytherea_center_points = easy_center_lines(cytherea_3d,20,sampling_freq,points_here_thresh);
caroliana_center_points = easy_center_lines(caroliana_3d,20,sampling_freq,points_here_thresh);

madracis6m_center_points = easy_center_lines(madracis6m_3d,15,5,points_here_thresh);
madracis15m_center_points = easy_center_lines(madracis15m_3d,15,5,points_here_thresh);
madracis20m_center_points = easy_center_lines(madracis20m_3d,15,5,points_here_thresh);

%% get rid of the bad center lines
millepora_center_points = cull_centers2(millepora_center_points,15);
loripes_center_points = cull_centers2(loripes_center_points,15);
cytherea_center_points = cull_centers2(cytherea_center_points,15);
caroliana_center_points = cull_centers2(caroliana_center_points,15);

madracis6m_center_points = cull_centers2(madracis6m_center_points,10);
madracis15m_center_points = cull_centers2(madracis15m_center_points,10);
madracis20m_center_points = cull_centers2(madracis20m_center_points,10);


%% get statistics from center lines
sampling_freq = 1;
thickness_samp = 5;

millepora_center_stats = center_line_analysis(millepora_3d,millepora_center_points,sampling_freq,thickness_samp);
loripes_center_stats = center_line_analysis(loripes_3d,loripes_center_points,sampling_freq,thickness_samp);
cytherea_center_stats = center_line_analysis(cytherea_3d,cytherea_center_points,sampling_freq,thickness_samp);
caroliana_center_stats = center_line_analysis(caroliana_3d,caroliana_center_points,sampling_freq,thickness_samp);

madracis6m_center_stats = center_line_analysis(madracis6m_3d,madracis6m_center_points,sampling_freq,thickness_samp);
madracis15m_center_stats = center_line_analysis(madracis15m_3d,madracis15m_center_points,sampling_freq,thickness_samp);
madracis20m_center_stats = center_line_analysis(madracis20m_3d,madracis20m_center_points,sampling_freq,thickness_samp);
%% and measure the branch angles

[millepora_brangles,millepora_brlengths] = spline_branch_angles2(millepora_branching_points_3d,millepora_branched_flags,millepora_center_stats,10);

[loripes_brangles,loripes_brlengths] = spline_branch_angles2(loripes_branching_points_3d,loripes_branched_flags,loripes_center_stats,10);

[cytherea_brangles,cytherea_brlengths] = spline_branch_angles2(cytherea_branching_points_3d,cytherea_branched_flags,cytherea_center_stats,10);

[caroliana_brangles,caroliana_brlengths] = spline_branch_angles2(caroliana_branching_points_3d,caroliana_branched_flags,caroliana_center_stats,10);

[madracis6m_brangles,madracis6m_brlengths] = spline_branch_angles2(madracis6m_branching_points_3d,madracis6m_branched_flags,madracis6m_center_stats,10);

[madracis15m_brangles,madracis15m_brlengths] = spline_branch_angles2(madracis15m_branching_points_3d,madracis15m_branched_flags,madracis15m_center_stats,10);

[madracis20m_brangles,madracis20m_brlengths] = spline_branch_angles2(madracis20m_branching_points_3d,madracis20m_branched_flags,madracis20m_center_stats,10);

%% surface areas, volumes, footprints, etc all in cm
[millepora_surface_area,millepora_volume] = sa_and_vol(millepora_downsampled, scale_ratio, millepora_scale/1e4);
[loripes_surface_area,loripes_volume] = sa_and_vol(loripes_downsampled, scale_ratio, loripes_scale/1e4);
[cytherea_surface_area,cytherea_volume] = sa_and_vol(cytherea_downsampled, scale_ratio, cytherea_scale/1e4);
[caroliana_surface_area,caroliana_volume] = sa_and_vol(caroliana_downsampled, scale_ratio, caroliana_scale/1e4);

[madracis6m_surface_area,madracis6m_volume] = sa_and_vol(madracis6m_downsampled, scale_ratio, madracis6m_scale/1e4);
[madracis15m_surface_area,madracis15m_volume] = sa_and_vol(madracis15m_downsampled, scale_ratio, madracis15m_scale/1e4);
[madracis20m_surface_area,madracis20m_volume] = sa_and_vol(madracis20m_downsampled, scale_ratio, madracis20m_scale/1e4);


[millepora_footprint_points, millepora_footprint_area] = get_footprint(millepora_3d, 3, millepora_scale/1e4);
[loripes_footprint_points, loripes_footprint_area] = get_footprint(loripes_3d, 3, loripes_scale/1e4);
[cytherea_footprint_points, cytherea_footprint_area] = get_footprint(cytherea_3d, 3, cytherea_scale/1e4);
[caroliana_footprint_points, caroliana_footprint_area] = get_footprint(caroliana_3d, 3, caroliana_scale/1e4);

[madracis6m_footprint_points, madracis6m_footprint_area] = get_footprint(madracis6m_3d, 3, madracis6m_scale/1e4);
[madracis15m_footprint_points, madracis15m_footprint_area] = get_footprint(madracis15m_3d, 3, madracis15m_scale/1e4);
[madracis20m_footprint_points, madracis20m_footprint_area] = get_footprint(madracis20m_3d, 3, madracis20m_scale/1e4);

[millepora_convhull_points, millepora_enclosing_volume] = get_enclosing_volume(millepora_3d, millepora_scale/1e4);
[loripes_convhull_points, loripes_enclosing_volume] = get_enclosing_volume(loripes_3d, loripes_scale/1e4);
[cytherea_convhull_points, cytherea_enclosing_volume] = get_enclosing_volume(cytherea_3d, cytherea_scale/1e4);
[caroliana_convhull_points, caroliana_enclosing_volume] = get_enclosing_volume(caroliana_3d, caroliana_scale/1e4);

[madracis6m_convhull_points, madracis6m_enclosing_volume] = get_enclosing_volume(madracis6m_3d, madracis6m_scale/1e4);
[madracis15m_convhull_points, madracis15m_enclosing_volume] = get_enclosing_volume(madracis15m_3d, madracis15m_scale/1e4);
[madracis20m_convhull_points, madracis20m_enclosing_volume] = get_enclosing_volume(madracis20m_3d, madracis20m_scale/1e4);


%%
diff_thresh = 5;

[millepora_deriv_means,millepora_deriv_variances,millepora_thicks_encountered,millepora_nn_dists] = centers_plane_pass(millepora_center_stats,millepora_3d, diff_thresh);
[loripes_deriv_means,loripes_deriv_variances,loripes_thicks_encountered,loripes_nn_dists] = centers_plane_pass(loripes_center_stats,loripes_3d, diff_thresh);
[cytherea_deriv_means,cytherea_deriv_variances,cytherea_thicks_encountered,cytherea_nn_dists] = centers_plane_pass(cytherea_center_stats,cytherea_3d, diff_thresh);
[caroliana_deriv_means,caroliana_deriv_variances,caroliana_thicks_encountered,caroliana_nn_dists] = centers_plane_pass(caroliana_center_stats,caroliana_3d, diff_thresh);

[madracis6m_deriv_means,madracis6m_deriv_variances,madracis6m_thicks_encountered,madracis6m_nn_dists] = centers_plane_pass(madracis6m_center_stats,madracis6m_3d, diff_thresh);
[madracis15m_deriv_means,madracis15m_deriv_variances,madracis15m_thicks_encountered,madracis15m_nn_dists] = centers_plane_pass(madracis15m_center_stats,madracis15m_3d, diff_thresh);
[madracis20m_deriv_means,madracis20m_deriv_variances,madracis20m_thicks_encountered,madracis20m_nn_dists] = centers_plane_pass(madracis20m_center_stats,madracis20m_3d, diff_thresh);

%% standardized surface areas per enclosing volume
n_cubes = 100;
cube_size = 40000;
tic
[millepora_surf_areas] = std_sa_encvol(millepora_downsampled,scale_ratio,millepora_scale,n_cubes,cube_size);
[loripes_surf_areas] = std_sa_encvol(loripes_downsampled,scale_ratio,loripes_scale,n_cubes,cube_size);
[cytherea_surf_areas] = std_sa_encvol(cytherea_downsampled,scale_ratio,cytherea_scale,n_cubes,cube_size);
[caroliana_surf_areas] = std_sa_encvol(caroliana_downsampled,scale_ratio,caroliana_scale,n_cubes,cube_size);
toc
tic
[madracis6m_surf_areas] = std_sa_encvol(madracis6m_downsampled,scale_ratio,madracis6m_scale,n_cubes,cube_size);
[madracis15m_surf_areas] = std_sa_encvol(madracis15m_downsampled,scale_ratio,madracis15m_scale,n_cubes,cube_size);
[madracis20m_surf_areas] = std_sa_encvol(madracis20m_downsampled,scale_ratio,madracis20m_scale,n_cubes,cube_size);
toc
tic
[sm_surf_areas] = std_sa_encvol(sm_outer_dense_slices,sm_scale_ratio,sm_um_pixel,n_cubes,cube_size,sm_branched_flags);
[cc297_surf_areas] = std_sa_encvol(cc297_outer_dense_slices,cc297_scale_ratio,cc297_um_pixel,n_cubes,cube_size,cc297_branched_flags);
[labrador_surf_areas] = std_sa_encvol(labrador_outer_dense_slices,labrador_scale_ratio,labrador_um_pixel,n_cubes,cube_size,labrador_branched_flags);
toc
%%
all_surf_data = [millepora_surf_areas';loripes_surf_areas';cytherea_surf_areas';caroliana_surf_areas';...
    madracis6m_surf_areas';madracis15m_surf_areas';madracis20m_surf_areas';sm_surf_areas';cc297_surf_areas';labrador_surf_areas';lds_surf_areas']/cube_size^3;
labels = [ones(n_cubes,1);...
    2*ones(n_cubes,1);...
    3*ones(n_cubes,1);...
    5*ones(n_cubes,1);...
    6*ones(n_cubes,1);...
    7*ones(n_cubes,1);...
    8*ones(n_cubes,1);...
    9*ones(n_cubes,1);...
    10*ones(n_cubes,1);...
    11*ones(n_cubes,1);...
    13*ones(n_cubes,1)];
    

boxplot(all_surf_data, labels,'orientation','horizontal','symbol','')
xlabel('Surface area / enclosing volume')


%%
figure();
subplot(2,2,1)
branch_angles = caroliana_brangles(:,:,4);
histogram(unique(branch_angles(branch_angles~=0 & branch_angles<90)),10)
title('A. caroliana')
xlabel('branch angle')
subplot(2,2,2)
branch_angles = cytherea_brangles(:,:,4);
histogram(unique(branch_angles(branch_angles~=0 & branch_angles<90)),10)
title('A. cytherea')
xlabel('branch angle')
subplot(2,2,3)
branch_angles = loripes_brangles(:,:,4);
histogram(unique(branch_angles(branch_angles~=0 & branch_angles<90)),10)
title('A. loripes')
xlabel('branch angle')
subplot(2,2,4)
branch_angles = millepora_brangles(:,:,4);
histogram(unique(branch_angles(branch_angles~=0 & branch_angles<90)),10)
title('A. millepora')
xlabel('branch angle')

figure();
subplot(3,1,1)
branch_angles = madracis6m_brangles(:,:,4);
histogram(unique(branch_angles(branch_angles~=0 & branch_angles<90)),10)
title('Madracis 6m')
xlabel('branch angle')
xlim([0,90])
subplot(3,1,2)
branch_angles = madracis15m_brangles(:,:,4);
histogram(unique(branch_angles(branch_angles~=0 & branch_angles<90)),10)
title('Madracis 15m')
xlabel('branch angle')
xlim([0,90])
subplot(3,1,3)
branch_angles = madracis20m_brangles(:,:,4);
histogram(unique(branch_angles(branch_angles~=0 & branch_angles<90)),10)
title('Madracis 20m')
xlabel('branch angle')
xlim([0,90])
