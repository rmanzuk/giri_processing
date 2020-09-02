function[inner_center_points,outer_center_points] = svd_center_lines(inner_3d,outer_3d,sampling_resolution,plt)
% This function takes 3d archaeo point clouds (densified or not) finds the
% center line of each archaeo at the disired sampling resolution. To do so,
% it tracks the axis of maximal variance (from SVD) and takes the centroid
% of the point cloud in a given sampling interval along that axis.
%
% IN
% inner_3d:: 1xn_archaeos cell array containing the densified or non-densified 3d outputs for
% inner clicked data from the densify_3d or make_clicking_3d functions.
%
% outer_3: 1xn_archaeos cell array containing the densified or non-densified 3d outputs for
% outer clicked data from the densify_3d or make_clicking_3d functions.
%
% sampling_resolution: interval (in the units of the coordinate space) at
% which you would like to sample for the center line points. So if your
% coordinate space is in microns, and you want to sample every 5 microns,
% this variable should be 5.
%
% plt: logical flag if the user would like the 3D center lines plotted. 
% 1 for plot, 0 for don't plot
% 
% OUT
% inner_center_points: 1xn_archaeo cell array where each cell contains the 3D
% coordinates for the center line of the inner tube for the given archaeo.
%
% outer_center_points: 1xn_archaeo cell array where each cell contains the 3D
% coordinates for the center line of the inner tube for the given archaeo.
%
% Ryan A. Manzuk 09/01/2020

    %% begin the function
    inner_center_points = {};
    for i = 1:numel(inner_3d)
            % extract 3d (dense) data for this archaeo
            current_points = inner_3d{i};
            % set up xyz data matrix
            xyz = current_points(:,1:3);
            % center the data
            xyz_mean = mean(xyz);
            xyz_centered = bsxfun(@minus, xyz, xyz_mean);
            % take the svd to understand how to rotate the archaeo
            [U,S,V] = svd(xyz_centered);
            % and rotate the data with respect to direction of maximal variance
            xyz_rotated = xyz_centered * V;
            % we'll need the range of rotated data along the primary access
            prim_range = range(xyz_rotated(:,1));
            % based upon range and min/max of primary axis, and sampling
            % resolution...where do we sample
            sample_ranges = linspace(min(xyz_rotated(:,1)),max(xyz_rotated(:,1)),round(prim_range/sampling_resolution)+2);
            %now go through and gather the centroids of all the ranges
            center_points = [];
            for j = 1:length(sample_ranges)-1
                in_range_logical = xyz_rotated(:,1) < sample_ranges(j+1) & xyz_rotated(:,1) >=sample_ranges(j);
                points_here = xyz_rotated(in_range_logical,:);
                center_points(j,:) = mean(points_here);
            end
            inner_center_points{i} = bsxfun(@plus, center_points*inv(V), xyz_mean); center_points*inv(V);
    end

    outer_center_points = {};
    for i = 1:numel(outer_3d)
            % extract 3d (dense) data for this archaeo
            current_points = outer_3d{i};
            % set up xyz data matrix
            xyz = current_points(:,1:3);
            % center the data
            xyz_mean = mean(xyz);
            xyz_centered = bsxfun(@minus, xyz, xyz_mean);
            % take the svd to understand how to rotate the archaeo
            [~,~,V] = svd(xyz_centered);
            % and rotate the data with respect to direction of maximal variance
            xyz_rotated = xyz_centered * V;
            % we'll need the range of rotated data along the primary access
            prim_range = range(xyz_rotated(:,1));
            % based upon range and min/max of primary axis, and sampling
            % resolution...where do we sample
            sample_ranges = linspace(min(xyz_rotated(:,1)),max(xyz_rotated(:,1)),round(prim_range/sampling_resolution)+2);
            %now go through and gather the centroids of all the ranges
            center_points = [];
            for j = 1:length(sample_ranges)-1
                in_range_logical = xyz_rotated(:,1) < sample_ranges(j+1) & xyz_rotated(:,1) >=sample_ranges(j);
                points_here = xyz_rotated(in_range_logical,:);
                center_points(j,:) = mean(points_here);
            end
            outer_center_points{i} = bsxfun(@plus, center_points*inv(V), xyz_mean); center_points*inv(V);
    end


    if plt
        for i = 1:numel(inner_center_points)
            plot3(inner_center_points{i}(:,1),inner_center_points{i}(:,2),inner_center_points{i}(:,3),'r')
            hold on
            plot3(outer_center_points{i}(:,1),outer_center_points{i}(:,2),outer_center_points{i}(:,3),'b')
        end
    else
        % do nothing
    end

end