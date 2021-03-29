function[center_points] = easy_center_lines(outlines,sampling_resolution,sampling_freq,points_here_thresh)
% This function takes 3d archaeo point clouds (densified or not) finds the
% center line of each archaeo at the disired sampling resolution, with the
% point clouds as is with no rotations.
%
% IN
% outlines: 1xn_branch cell array containing the densified or non-densified
% 3d outputs for outline points from the densify_3d or make_clicking_3d
% functions.
%
% sampling_resolution: interval (in the units of the coordinate space)
% above and below your sampling point over which you would like to consider
% for the center calculation. For example, if you would like to consider
% points 5 units above and below each sampled point when determining the
% center, this value should be 5.
%
% sampling freq: interval (in the units of the coordinate space) at
% which you would like to sample for the center line points. So if your
% coordinate space is in microns, and you want to sample every 5 microns,
% this variable should be 5.
%
% points_here_thresh: minimum number of points that need to be within an
% interval for the function to assess the centroid for the interval.
%
% 
% OUT
% center_points: 1xn_branch cell array where each cell contains the 3D
% coordinates for the center line of the given branch.
%
% Ryan A. Manzuk 09/08/2020

    %% begin the function
    center_points = {};
    for i = 1:numel(outlines)
            % extract 3d (dense) data for this archaeo
            current_points = outlines{i};
            % set up xyz data matrix
            xyz = current_points(:,1:3);
            % we'll need the range of 3d data along the z axis
            prim_range = range(xyz(:,3));
            % based upon range and min/max of primary axis, and sampling
            % frequency...where do we sample
            sampled_points = linspace(min(xyz(:,3)),max(xyz(:,3)),round(prim_range/sampling_freq)+2);
            %now go through and gather the centroids of all the ranges
            branch_centers = [];
            for j = 1:length(sampled_points)
                % which points are in the range of this sample point +/-
                % the sampling resolution
                in_range_logical = xyz(:,3) < sampled_points(j)+sampling_resolution & xyz(:,3) >sampled_points(j)-sampling_resolution;
                if sum(in_range_logical) >= points_here_thresh
                    points_here = xyz(in_range_logical,:);
                    branch_centers(j,:) = [mean(points_here(:,1:2)),sampled_points(j)];
                else
                    branch_centers(j,:) = [NaN,NaN,NaN];
                end
            end
            center_points{i} = branch_centers;
    end
end