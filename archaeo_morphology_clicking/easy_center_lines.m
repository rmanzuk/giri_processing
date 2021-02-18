function[inner_center_points,outer_center_points] = easy_center_lines(inner_3d,outer_3d,sampling_resolution,sampling_freq,points_here_thresh)
% This function takes 3d archaeo point clouds (densified or not) finds the
% center line of each archaeo at the disired sampling resolution, with the
% point clouds as is with no rotations.
%
% IN
% inner_3d:: 1xn_archaeos cell array containing the densified or non-densified 3d outputs for
% inner clicked data from the densify_3d or make_clicking_3d functions.
%
% outer_3d: 1xn_archaeos cell array containing the densified or non-densified 3d outputs for
% outer clicked data from the densify_3d or make_clicking_3d functions.
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
% inner_center_points: 1xn_archaeo cell array where each cell contains the 3D
% coordinates for the center line of the inner tube for the given archaeo.
%
% outer_center_points: 1xn_archaeo cell array where each cell contains the 3D
% coordinates for the center line of the inner tube for the given archaeo.
%
% Ryan A. Manzuk 09/08/2020

    %% begin the function
    inner_center_points = {};
    for i = 1:numel(inner_3d)
            % extract 3d (dense) data for this archaeo
            current_points = inner_3d{i};
            % set up xyz data matrix
            xyz = current_points(:,1:3);
            % we'll need the range of 3d data along the z axis
            prim_range = range(xyz(:,3));
            % based upon range and min/max of primary axis, and sampling
            % frequency...where do we sample
            sampled_points = linspace(min(xyz(:,3)),max(xyz(:,3)),round(prim_range/sampling_freq)+2);
            %now go through and gather the centroids of all the ranges
            center_points = [];
            for j = 1:length(sampled_points)
                % which points are in the range of this sample point +/-
                % the sampling resolution
                in_range_logical = xyz(:,3) < sampled_points(j)+sampling_resolution & xyz(:,3) >sampled_points(j)-sampling_resolution;
                if sum(in_range_logical) >= points_here_thresh
                    points_here = xyz(in_range_logical,:);
                    center_points(j,:) = [mean(points_here(:,1:2)),sampled_points(j)];
                else
                    center_points(j,:) = [NaN,NaN,NaN];
                end
            end
            inner_center_points{i} = center_points;
    end

    outer_center_points = {};
    for i = 1:numel(outer_3d)
            % extract 3d (dense) data for this archaeo
            current_points = outer_3d{i};
            % set up xyz data matrix
            xyz = current_points(:,1:3);
            % we'll need the range of 3d data along the z axis
            prim_range = range(xyz(:,3));
            % based upon range and min/max of primary axis, and sampling
            % frequency...where do we sample
            sampled_points = linspace(min(xyz(:,3)),max(xyz(:,3)),round(prim_range/sampling_freq)+2);
            %now go through and gather the centroids of all the ranges
            center_points = [];
            for j = 1:length(sampled_points)
                % which points are in the range of this sample point +/-
                % the sampling resolution
                in_range_logical = xyz(:,3) < sampled_points(j)+sampling_resolution & xyz(:,3) >sampled_points(j)-sampling_resolution;
                if sum(in_range_logical) >= points_here_thresh
                    points_here = xyz(in_range_logical,:);
                    center_points(j,:) = [mean(points_here(:,1:2)),sampled_points(j)];
                else
                    center_points(j,:) = [NaN,NaN,NaN];
                end
            end
            outer_center_points{i} = center_points;
    end

end