function[center_points] = svd_center_lines(outlines_3d,sampling_resolution,points_here_thresh)
% This function takes 3d archaeo point clouds (densified or not) finds the
% center line of each archaeo at the disired sampling resolution. To do so,
% it tracks the axis of maximal variance (from SVD) and takes the centroid
% of the point cloud in a given sampling interval along that axis.
%
% IN
% outlines_3d:: 1xn_branch cell array containing the densified or non-densified 3d outputs for
% clicked data from the densify_3d or make_clicking_3d functions.
%
%
% sampling_resolution: interval (in the units of the coordinate space) at
% which you would like to sample for the center line points. So if your
% coordinate space is in microns, and you want to sample every 5 microns,
% this variable should be 5.
%
% points_here_thresh: minimum number of points that need to be within an
% interval for the function to assess the centroid for the interval.
%
% plt: logical flag if the user would like the 3D center lines plotted. 
% 1 for plot, 0 for don't plot
% 
% OUT
% center_points: 1xn_branch cell array where each cell contains the 3D
% coordinates for the center line of the inner tube for the given branch.
%
% Ryan A. Manzuk 09/01/2020

    %% begin the function
    center_points = {};
    for i = 1:numel(outlines_3d)
        % extract 3d (dense) data for this archaeo
        current_points = outlines_3d{i};
        % set up xyz data matrix
        xyz = current_points(:,1:3);
        % center the data
        xyz_mean = mean(xyz);
        xyz_centered = bsxfun(@minus, xyz, xyz_mean);
        % take the svd to understand how to rotate the archaeo
        [~,~,V] = svd(xyz_centered,'econ');
        % and rotate the data with respect to direction of maximal variance
        xyz_rotated = xyz_centered * V;
        % we'll need the range of rotated data along the primary access
        prim_range = range(xyz_rotated(:,1));
        % based upon range and min/max of primary axis, and sampling
        % resolution...where do we sample
        sample_ranges = linspace(min(xyz_rotated(:,1)),max(xyz_rotated(:,1)),round(prim_range/sampling_resolution)+2);
        %now go through and gather the centroids of all the ranges
        these_center_points = [];
        for j = 1:length(sample_ranges)-1
            in_range_logical = xyz_rotated(:,1) < sample_ranges(j+1) & xyz_rotated(:,1) >=sample_ranges(j);
            if sum(in_range_logical) >= points_here_thresh
                points_here = xyz_rotated(in_range_logical,:);
                these_center_points(j,:) = mean(points_here);
            else
                these_center_points(j,:) = [NaN,NaN,NaN];
            end
        end
        centers_moved_back = bsxfun(@plus, these_center_points*inv(V), xyz_mean); these_center_points*inv(V);
        if centers_moved_back(1,3) > centers_moved_back(end,3)
            centers_moved_back = flip(centers_moved_back);
        end
        center_points{i} = centers_moved_back;
    end
end