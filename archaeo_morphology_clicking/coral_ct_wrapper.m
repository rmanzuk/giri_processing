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
colors = [get(groot,'defaultAxesColorOrder');get(groot,'defaultAxesColorOrder');get(groot,'defaultAxesColorOrder')];

for i = 1:numel(final_outers)
    for j = 1:numel(final_outers{i})
        if ~isempty(final_outers{i}{j})
            scatter3(final_outers{i}{j}(:,1),final_outers{i}{j}(:,2),(-5*i).*ones(length(final_outers{i}{j}(:,2)),1),40,colors(j,:))
        else
        end
            hold on
    end
end