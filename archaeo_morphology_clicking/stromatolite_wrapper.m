% script to run through and make a few stromatolites, coalesc them and
% measure their params
%
% R. A. Manzuk 03/29/2021
%% set up the general stuff

% these params hold for all stromatolites
height = 10;
vert_sampling_rate = 10;
circle_sampling_rate = 10;
bump_size = 1;
stromat_scale_ratio = 1/vert_sampling_rate;

% unique set of radii for each stromatolite
rads1 = [5,6,7,8,9,10,10,10,10,10];
rads2 = [6,6.5,8,8,9,10,10,10,10,10];
rads3 = [2,4,6,8,9,10,10,10,9,8];
rads4 = [1,2,3,4,5,6,7,8,9,10];

%% Make the stromatolites
[stromat1] = make_stromatolite(rads1,height,vert_sampling_rate,circle_sampling_rate,bump_size);
[stromat2] = make_stromatolite(rads2,height,vert_sampling_rate,circle_sampling_rate,bump_size);
[stromat3] = make_stromatolite(rads3,height,vert_sampling_rate,circle_sampling_rate,bump_size);
[stromat4] = make_stromatolite(rads4,height,vert_sampling_rate,circle_sampling_rate,bump_size);

%% now translate the stromatolites so they can be in different positions and we can coalesc them
translation_length = 12;

% move stromat2 in positive x direction
for i = 1:numel(stromat2)
    points_here = stromat2{i};
    points_here(:,1) = points_here(:,1) + translation_length;
    stromat2{i} = points_here;
end

% move stromat3 in positive y direction
for i = 1:numel(stromat3)
    points_here = stromat3{i};
    points_here(:,2) = points_here(:,2) + translation_length;
    stromat3{i} = points_here;
end

% move stromat4 in positive x and y direction
for i = 1:numel(stromat4)
    points_here = stromat4{i};
    points_here(:,1) = points_here(:,1) + translation_length;
    points_here(:,2) = points_here(:,2) + translation_length;
    stromat4{i} = points_here;
end

%% coalesc
all_stromats = {stromat1,stromat2,stromat3,stromat4};

%% make 3d and measure

% identify branch points
[stromat_branched_flags,stromat_branch_points_3d] = process_branched_network(all_stromats,stromat_scale_ratio);

% make those stromats 3d
stromats_3d = make_clicking_3d(all_stromats,stromat_scale_ratio);

% take the center lines
sampling_resolution = 1;
sampling_freq = 1;
points_here_thresh = 20;
stromat_center_points = easy_center_lines(stromats_3d,sampling_resolution,sampling_freq,points_here_thresh);

% and get some data
thickness_samp = 0.5;
diff_thresh = 5;

stromat_center_stats = center_line_analysis(stromats_3d,stromat_center_points,sampling_freq,thickness_samp);

[stromat_surface_area,stromat_volume] = sa_and_vol(all_stromats,stromat_scale_ratio,1,stromat_branched_flags);
[stromat_footprint_points, stromat_footprint_area] = get_footprint(stromats_3d, 3, 1);
[stromat_convhull_points, stromat_enclosing_volume] = get_enclosing_volume(stromats_3d, 1);
