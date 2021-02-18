% scrip to hopefully automate the tracing of edges in coral ct scans
% (analagous to manual clicking in grinder stacks) to then put traces
% through same measurement pipeline.

% R. A. Manzuk 12/15/2020
%% set up the image directory

input_folder = '/Users/rmanzuk/Desktop/coral_ct_slices';

file_pattern = fullfile(input_folder, '*.tif');
tifs = dir(file_pattern);
base_names = natsortfiles({tifs.name});
fulln=numel(base_names);
%% loop though all images and extract edges

% tunable parameters
blur_kernel_size = 10;
island_size = 200;
hole_size = 200;
distance_thresh = 10;

% set up final_outers to receive edges
final_outers = {};
% set up a counter
counter = 0;
for i = 1:numel(base_names)
    %read the image
    this_im = imread(fullfile(input_folder, base_names{i}));
    % use function to get this slice's edges
    edge_coords = get_slice_edges(this_im,blur_kernel_size,island_size,hole_size);
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

load('/Users/rmanzuk/Desktop/branch_angle_project/coral_ct_data/millepora_cleaned.mat');
load('/Users/rmanzuk/Desktop/branch_angle_project/coral_ct_data/loripes_cleaned.mat');
load('/Users/rmanzuk/Desktop/branch_angle_project/coral_ct_data/cytherea_cleaned.mat');
load('/Users/rmanzuk/Desktop/branch_angle_project/coral_ct_data/caroliana_cleaned.mat');

%%
% rearrange the cell array to be 1xn_branches
millepora_reshaped = reshape_coral_cell(millepora_cleaned);
loripes_reshaped = reshape_coral_cell(loripes_cleaned);
cytherea_reshaped = reshape_coral_cell(cytherea_cleaned);
caroliana_reshaped = reshape_coral_cell(caroliana_cleaned);
clear millepora_cleaned
clear loripes_cleaned
clear cytherea_cleaned
clear caroliana_cleaned

% we also need the edges sorted as if going around the circle, not with
% sequential indices
millepora_resorted = sort_outline_points(millepora_reshaped);
loripes_resorted = sort_outline_points(loripes_reshaped);
cytherea_resorted = sort_outline_points(cytherea_reshaped);
caroliana_resorted = sort_outline_points(caroliana_reshaped);
clear millepora_reshaped
clear loripes_reshaped
clear cytherea_reshaped
clear caroliana_reshaped

% and then let's downsample the slices
densify_factor = 0.1;

[~,millepora_downsampled] = densify_slices(millepora_resorted,millepora_resorted,densify_factor);
[~,loripes_downsampled] = densify_slices(loripes_resorted,loripes_resorted,densify_factor);
[~,cytherea_downsampled] = densify_slices(cytherea_resorted,cytherea_resorted,densify_factor);
[~,caroliana_downsampled] = densify_slices(caroliana_resorted,caroliana_resorted,densify_factor);
clear millepora_resorted
clear loripes_resorted
clear cytherea_resorted
clear caroliana_resorted
%% The scales of each specimen
% caroliniana 
caroliana_scale = 81.1170;
% cytherea
cytherea_scale = 107.8580;
% loripes 
loripes_scale = 96.5800;
% millepora 
millepora_scale = 49.9990;
%% and get the branching points
distance_threshold = 20;
combining_threshold = 20;
scale_ratio = 1;
[millepora_branched_flags,millepora_branching_points_3d,millepora_combined] = id_branch_points(millepora_downsampled,scale_ratio,distance_threshold,combining_threshold);
[loripes_branched_flags,loripes_branching_points_3d,loripes_combined] = id_branch_points(loripes_downsampled,scale_ratio,distance_threshold,combining_threshold);
[cytherea_branched_flags,cytherea_branching_points_3d,cytherea_combined] = id_branch_points(cytherea_downsampled,scale_ratio,distance_threshold,combining_threshold);
[caroliana_branched_flags,caroliana_branching_points_3d,caroliana_combined] = id_branch_points(caroliana_downsampled,scale_ratio,distance_threshold,combining_threshold);
clear millepora_downsampled
clear loripes_downsampled
clear cytherea_downsampled
clear caroliana_downsampled
%% take slice data and make 3d
[~,millepora_3d] = make_clicking_3d(millepora_combined,millepora_combined,scale_ratio,0);
[~,loripes_3d] = make_clicking_3d(loripes_combined,loripes_combined,scale_ratio,0);
[~,cytherea_3d] = make_clicking_3d(cytherea_combined,cytherea_combined,scale_ratio,0);
[~,caroliana_3d] = make_clicking_3d(caroliana_combined,caroliana_combined,scale_ratio,0);
clear millepora_combined
clear loripes_combined
clear cytherea_combined
clear caroliana_combined
%% get the center lines
sampling_resolution = 30;
sampling_freq = 1;
points_here_thresh = 20;
[~,millepora_center_points] = easy_center_lines(millepora_3d,millepora_3d,sampling_resolution,sampling_freq,points_here_thresh);
[~,loripes_center_points] = easy_center_lines(loripes_3d,loripes_3d,sampling_resolution,sampling_freq,points_here_thresh);
[~,cytherea_center_points] = easy_center_lines(cytherea_3d,cytherea_3d,sampling_resolution,sampling_freq,points_here_thresh);
[~,caroliana_center_points] = easy_center_lines(caroliana_3d,caroliana_3d,sampling_resolution,sampling_freq,points_here_thresh);
%% get statistics from center lines
sampling_freq = 1;
thickness_samp = 5;

[~,millepora_center_stats] = center_line_analysis(millepora_3d,millepora_3d,millepora_center_points,sampling_freq,thickness_samp);
[~,loripes_center_stats] = center_line_analysis(loripes_3d,loripes_3d,loripes_center_points,sampling_freq,thickness_samp);
[~,cytherea_center_stats] = center_line_analysis(cytherea_3d,cytherea_3d,cytherea_center_points,sampling_freq,thickness_samp);
[~,caroliana_center_stats] = center_line_analysis(caroliana_3d,caroliana_3d,caroliana_center_points,sampling_freq,thickness_samp);
%% and measure the branch angles
% lets do it considering multiple lenghts for the calculation
lengths_considered = [0.1,0.5,1,2,3]; % in centimeters
% and apply the scale of the samples
millepora_lengths = lengths_considered./(millepora_scale/1e4);
loripes_lengths = lengths_considered./(loripes_scale/1e4);
cytherea_lengths = lengths_considered./(cytherea_scale/1e4);
caroliana_lengths = lengths_considered./(caroliana_scale/1e4);

% analyze each sample at each length
millepora_brangles = [];
millepora_radii = [];
for i = 1: length(millepora_lengths)
    [millepora_brangles(:,:,i),millepora_radii(:,:,i)] = spline_branch_angles(millepora_branching_points_3d,millepora_branched_flags,millepora_center_stats,millepora_lengths(i));
end

loripes_brangles = [];
loripes_radii = [];
for i = 1: length(loripes_lengths)
    [loripes_brangles(:,:,i),loripes_radii(:,:,i)] = spline_branch_angles(loripes_branching_points_3d,loripes_branched_flags,loripes_center_stats,loripes_lengths(i));
end

cytherea_brangles = [];
cytherea_radii = [];
for i = 1: length(cytherea_lengths)
    [cytherea_brangles(:,:,i),cytherea_radii(:,:,i)] = spline_branch_angles(cytherea_branching_points_3d,cytherea_branched_flags,cytherea_center_stats,cytherea_lengths(i));
end

caroliana_brangles = [];
caroliana_radii = [];
for i = 1: length(caroliana_lengths)
    [caroliana_brangles(:,:,i),caroliana_radii(:,:,i)] = spline_branch_angles(caroliana_branching_points_3d,caroliana_branched_flags,caroliana_center_stats,caroliana_lengths(i));
end

%%
figure();
for j=1:length(lengths_considered);
    subplot(3,2,j)
    branch_angles = caroliana_brangles(:,:,j);
    histogram(unique(branch_angles(branch_angles~=0 & branch_angles<90)),10)
    title([num2str(lengths_considered(j)) ' cm']);
end