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
%%
% Cleaning up the automated tracing output

% first, get rid of edges that don't span more than ~5 slices
slice_cutoff = 3;

cleaned_outers = remove_spurious_edges(final_outers, slice_cutoff);

% and we also need to rearrange the cell array to be 1xn_branches
reshaped_outers = reshape_coral_cell(cleaned_outers);

% we also need the edges sorted as if going around the circle, not with
% sequential indices
resorted_outers = sort_outline_points(reshaped_outers);

% and then let's downsample the slices
densify_factor = 0.1;
[~,downsampled_outers] = densify_slices(resorted_outers,resorted_outers,densify_factor);

%% take slice data and make 3d
scale_ratio = 10;
[~,outers_3d] = make_clicking_3d(downsampled_outers,downsampled_outers,scale_ratio,1);
%% and get the branching points
distance_threshold = 50;
[branched_flags,branching_points_3d] = id_branch_points(downsampled_outers,scale_ratio,distance_threshold);
%% get the center lines
sampling_resolution = 20;
points_here_thresh = 10;
[~,center_points] = easy_center_lines(outers_3d,outers_3d,sampling_resolution,points_here_thresh);
%% get statistics from center lines
sampling_freq = 1;
thickness_samp = 3;
[~,outer_center_stats] = center_line_analysis(outers_3d,outers_3d,center_points,sampling_freq,thickness_samp);
%% and measure the branch angles
branch_angles = spline_branch_angles(branching_points_3d,outer_center_stats,300);