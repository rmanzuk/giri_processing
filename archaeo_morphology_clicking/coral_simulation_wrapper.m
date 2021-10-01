% wrapper script to make some measurements on all of the simulated growth
% experiments
%
% R. A. Manzuk 07/13/2021
%%
% we'll just loop through all of the stacks to save some time
% so set that up
main_directory = '/Users/ryan/Desktop/branch_angle_project/coral_ct_data/Simulations/tif_stacks';
subdirectories = dir(fullfile(main_directory,'*'));
% list of subfolders of main_directory
subs_list = setdiff({subdirectories([subdirectories.isdir]).name},{'.','..'});

% set up cells for final measurements
simulation_brangles = {};
simulation_brlengths = {};
thicks_encountered = {};
nn_dists = {};
deriv_means = {};
deriv_variances = {};
for i = 1:numel(subs_list)
    tic
    % id the stack folder
    file_pattern = fullfile([main_directory '/' subs_list{i}], '*.tiff');
    this_stack = dir(file_pattern);
    base_names = natsortfiles({this_stack.name});
    % all stacks are upside down, so read backwards
    base_names = flip(base_names);
    % set up cell to get outlines and a counter
    final_outers = {};
    counter = 0;
    final_binaries = [];
    % iterate through all images to get the outlines
    for j = 1:numel(base_names)
        this_im = imread(fullfile([main_directory '/' subs_list{i}], base_names{j}));
        % make it a double
        im_double = im2double(this_im);
        % get the slice edges, we expect islands to be pretty small but
        % holes to be pretty big
        blur_kernel_size = 0.01;
        island_size = 100;
        hole_size = 500000;
        [edge_coords,final_binaries(:,:,j)] = get_slice_edges2(im_double,blur_kernel_size,island_size,hole_size);
        
        if ~isempty(edge_coords)
        % we have edges here, so update the counter
        counter = counter + 1;
            if counter == 1
                % if these are our first edges, we just take them, no questions
                % asked
                final_outers{j} = edge_coords;
            else
                % otherwise we have to match the edges, use the function
                distance_thresh = 50;
                final_outers{j} = match_edges(edge_coords,final_outers{j-1},distance_thresh);
            end
        else
            if counter == 0
                % if j is 1 and we dont have any edges we don't do
                % anything
            else
                % we need to just put empty cells for all of the old
                % branches
                counter = counter + 1;
                final_outers{j} = cell(1,numel(final_outers{j-1}));
            end
        end
    end
    % remove empty slices and reshape the outline array to be 1 x n_branches
    final_outers = final_outers(~cellfun('isempty',final_outers));
    outlines_reshaped = reshape_coral_cell(final_outers);
    % we also need the edges sorted as if going around the circle, not with
    % sequential indices
    outlines_resorted = sort_outline_points(outlines_reshaped);
    % and then let's downsample the slices
    densify_factor = 0.5;
 
    outlines_downsampled = densify_slices(outlines_resorted,densify_factor);
    % and get the branching points
    distance_threshold = 50;
    combining_threshold = 50;
    scale_ratio = 10;
    [branched_flags,branching_points_3d,outlines_combined] = id_branch_points(outlines_downsampled,scale_ratio,distance_threshold,combining_threshold);
    
    % take slice data and make 3d
    outlines_3d = make_clicking_3d(outlines_combined,scale_ratio);
    
    % get the center lines
    sampling_resolution = 50;
    sampling_freq = 0.5;
    points_here_thresh = 30;

    center_points = easy_center_lines(outlines_3d,sampling_resolution,sampling_freq,points_here_thresh);
    
    % get rid of the bad center lines
    length_thresh = 10;
    center_points = cull_centers2(center_points,length_thresh);
    
    % analyze the center lines
    sampling_freq = 1;
    thickness_samp = 5;

    center_stats = center_line_analysis(outlines_3d,center_points,sampling_freq,thickness_samp);
    
    % and measure the branch angles
    [simulation_brangles{i},simulation_brlengths{i}] = spline_branch_angles2(branching_points_3d,branched_flags,center_stats,10);
    % and get center stats
    [deriv_means{i},deriv_variances{i},thicks_encountered{i},nn_dists{i}] = centers_plane_pass(center_stats,outlines_3d, 11);
    toc
end

%% get out some statistics
brangle_means = [];
brangle_stdevs = [];
brangle_medians = [];
nn_means = [];
nn_stddevs = [];
thickness_means = [];
thickness_stddevs = [];
for i = 1:numel(simulation_brangles)
    these_brangles = simulation_brangles{i};
    these_nns = nn_dists{i};
    these_thicks = thicks_encountered{i};
    brangle_means(i) = mean(these_brangles(these_brangles > 0),'omitnan');
    brangle_stdevs(i) = std(these_brangles(these_brangles > 0),'omitnan');
    brangle_medians(i) = median(these_brangles(these_brangles > 0),'omitnan');
    nn_means(i) = mean(these_nns,'all','omitnan');
    nn_stdevs(i) = std(these_nns,[],'all','omitnan');
    thickness_means(i) = mean(these_thicks,'all','omitnan');
    thickness_stddevs(i) = std(these_thicks,[],'all','omitnan');
end
%% box plots of some stats for changing diffusive characteristics (ext. fig 6)

figure()
% assemble branching angles for box plot
all_brangle_data = [simulation_brangles{1}(simulation_brangles{1} >0);simulation_brangles{2}(simulation_brangles{2} >0)];
labels = [ones(numel(simulation_brangles{1}(simulation_brangles{1} >0)),1);...
    2*ones(numel(simulation_brangles{2}(simulation_brangles{2} >0)),1)];
    
subplot(1,2,1)
boxplot(all_brangle_data, labels,'symbol','')
ylabel('Branching Angle [degrees]')
ylim([0,90])

% assemble branch thicknesses for box plot
all_thickness_data = [thicks_encountered{1}(:);thicks_encountered{2}(:)];
labels = [ones(numel(thicks_encountered{1}(:)),1);...
    2*ones(numel(thicks_encountered{2}(:)),1)];
    
subplot(1,2,2)
boxplot(all_thickness_data, labels,'symbol','')
ylabel('Branch thickness [arbitrary units]')
ylim([0,100])

%% box plots of some stats when adding light sensitivity (ext. fig. 6)
 
figure()
% assemble branching angles for box plot
all_brangle_data = [simulation_brangles{1}(simulation_brangles{1} >0);simulation_brangles{12}(simulation_brangles{12} >0)];
labels = [ones(numel(simulation_brangles{1}(simulation_brangles{1} >0)),1);...
    2*ones(numel(simulation_brangles{12}(simulation_brangles{12} >0)),1)];
    
subplot(1,2,1)
boxplot(all_brangle_data, labels,'symbol','')
ylabel('Branching Angle [degrees]')
ylim([0,90])

% assemble branch headings
all_deriv_data = [deriv_variances{1}(:);deriv_variances{12}(:)];
labels = [ones(numel(deriv_variances{1}(:)),1);...
    2*ones(numel(deriv_variances{12}(:)),1)];
    
subplot(1,2,2)
boxplot(all_deriv_data, labels,'symbol','')
ylabel('Branch heading variance')
ylim([0,2])