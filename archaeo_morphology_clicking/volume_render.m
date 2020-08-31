function [inner_one_hot,outer_one_hot] = volume_render(inner_3d,outer_3d,mult_fac)
% This function takes 3d point clouds for inner and outer clicked archaeo
% data and turns them into volume renderings based upone which voxels are
% deemded to be inside the polygons defined by the point clouds.

% IN
% inner_3d: 1xn_archaeos cell array containing 3d point data for inner
% circles. Can be rotated, densified, whatever, just not set up slice-wise
%
% outer_3d:  1xn_archaeos cell array containing 3d point data for outer
% circles. Can be rotated, densified, whatever, just not set up slice-wise
%
% mult_fac: Factor by which the coordinate space should be multiplied to 
% increase volume rendering resolution....For example, if your point cloud
% occupies a space that is 200x300x400, a mult_fac of 2 will make the
% rendering matrix 400x600x800. Obviously difficult to have huge 3d
% matrices in memory, so hard to have a mult_fac greater than 2 for
% archaeos.
% 
% 
% OUT
% inner_one_hot: one-hot 3d matrix representing the volume rendering at 
% desired demensions of the inner 3d tubes. voxels with a 1 are part of the
% inner tube. Voxels with a 0 are not.
%
% outer_one_hot: one-hot 3d matrix representing the volume rendering at 
% desired demensions of the outer 3d tubes. voxels with a 1 are part of the
% outer tube. Voxels with a 0 are not.
%
% R. A. Manzuk, 08/31/2020
    %% begin the function
    all_inn = vertcat(inner_3d{[1:end]});
    all_out = vertcat(outer_3d{[1:end]});

    minima = min(vertcat(all_inn(:,1:3),all_out(:,1:3)));
    maxima = max(vertcat(all_inn(:,1:3),all_out(:,1:3)));

    ranges = maxima - minima;

    inner_one_hot = zeros(ceil(ranges(1))*mult_fac,ceil(ranges(2))*mult_fac,ceil(ranges(3))*mult_fac);
    outer_one_hot = zeros(ceil(ranges(1))*mult_fac,ceil(ranges(2))*mult_fac,ceil(ranges(3))*mult_fac);

    index_pairs = zeros(size(inner_one_hot,1)*size(inner_one_hot,2),2);
    for i = 1:size(inner_one_hot,1)
        for j = 1:size(inner_one_hot,2)
            index_pairs(((i-1)*size(inner_one_hot,1))+j,:) = [i,j];
        end
    end

    for i = 1:numel(inner_3d)
        current_points = inner_3d{i}(:,1:3);
        translated_points = (current_points-minima).*mult_fac;
        for z = floor(min(translated_points(:,3))):ceil(max(translated_points(:,3)))-1
            z_points = translated_points(:,3);
            logical_indices = z_points > z & z_points < z+1;
            slice_xy = translated_points(logical_indices,1:2);
            if size(slice_xy,1) >= 7
                conv_hull = convhull(slice_xy);
                in_shell = inpolygon(index_pairs(:,1),index_pairs(:,2),slice_xy(conv_hull,1),slice_xy(conv_hull,2));
                turn_on = [index_pairs(in_shell,:),ones(sum(in_shell),1).*z];
                for j = 1:sum(in_shell)
                    inner_one_hot(turn_on(j,1),turn_on(j,2),turn_on(j,3)) = 1;
                end
            else
            end
        end

    end

    for i = 1:numel(outer_3d)
        current_points = outer_3d{i}(:,1:3);
        translated_points = (current_points-minima).*mult_fac;
        for z = floor(min(translated_points(:,3))):ceil(max(translated_points(:,3)))-1
            z_points = translated_points(:,3);
            logical_indices = z_points > z & z_points < z+1;
            slice_xy = translated_points(logical_indices,1:2);
            if size(slice_xy,1) >= 5
                conv_hull = convhull(slice_xy);
                in_shell = inpolygon(index_pairs(:,1),index_pairs(:,2),slice_xy(conv_hull,1),slice_xy(conv_hull,2));
                turn_on = [index_pairs(in_shell,:),ones(sum(in_shell),1).*z];
                for j = 1:sum(in_shell)
                    outer_one_hot(turn_on(j,1),turn_on(j,2),turn_on(j,3)) = 1;
                end
            else
            end
        end

    end
end