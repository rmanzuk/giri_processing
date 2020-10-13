function [deriv_means, deriv_variances, thicks_encountered] = centers_plane_pass(outer_center_stats,outer_3d_rotated, diff_thresh)
% This function takes the derivative argument of the center stats output
% and conducts a rotation based upon the mean archaeo heading and the
% conducts a plane pass from bottom to top to understand the change in
% derivative and thereby things like parallelness and mean inclination.
%
% IN
% outer_center_stats: Structure containing 1xn_archaeo cell arrays where each 
% cell contains a given measure of the outer tube for the given archaeo.
% (output of the center_line_analysis.m function)
%
% outer_3d_rotated: 1xn_archaeo cell array where each cell contains the 3D
% coordinates for all outer clicked points (probably should be densified and 
% rotated for a given archaeo.
%
% diff_thresh: voxel unit within which points within an archaeo will be
% considdered to be part of the current plane. For example, if you would
% like each plane capture points +/- 5 voxels, this variable should be 5.
%
% OUT
% deriv_means: n_planes x 3 matrix giving the mean 3D derivative
% encountered by each plane going from bottom to top. n_planes is just
% dictated by the voxel span of the model in the vertical direction. If it
% covers 500 voxels, n_planes will be 500.
%
% deriv_variances: n_planes x 3 matrix giving the 3D variance of the derivative
% encountered by each plane going from bottom to top.
%
% thicks_encountered:  n_archaeos x n_planes matrix that can be analyzed
% for trends in mean thickness.
%
% R. A. Manzuk 10/13/2020
    %% begin the function
    % first, we need to extract all branch derivatives and find their mean
    all_derivs = [];
    for i = 1:numel(outer_center_stats.derivative)
        all_derivs = [all_derivs;outer_center_stats.derivative{i}];
    end
    mean_deriv = mean(all_derivs,1);
    
    % based on the mean derivative, find the inclination and declination
    % for rotation
    [declination,slope_run] = cart2pol(mean_deriv(1),mean_deriv(2));
    inclination = atand(mean_deriv(3)/slope_run);
    declination = rad2deg(declination);
    declination = -1.*(declination - 90);
    
    %fix declination if necessary
    if declination < 0
        declination = 360 + declination;
    end
    
    % set up rotation matrices
    z_rot_mat = [cosd(declination), -sind(declination), 0;
                sind(declination), cosd(declination), 0; 0, 0, 1];
    x_rot_mat = [1,0,0;0,cosd(90-inclination), -sind(90-inclination);
                0, sind(90-inclination), cosd(90-inclination)];
    
    % rotate everybody and centers based upon mean derivative
    outer_rotated_deriv = {};
    centers_rotated_deriv = {};
    for j = 1:numel(outer_3d_rotated)
        these_outers = outer_3d_rotated{j}(:,1:3);
        these_centers = outer_center_stats.spline{j};

        rot1_outers = z_rot_mat * these_outers';
        rot1_centers = z_rot_mat * these_centers';

        outer_rotated_deriv{j} = (x_rot_mat * rot1_outers)';
        centers_rotated_deriv{j} = (x_rot_mat * rot1_centers)';
    end
    
    % we need to know the range of z values covered by the newly rotated
    % model
    all_z = [];
    for k = 1:numel(outer_rotated_deriv)
        all_z = [all_z;outer_rotated_deriv{k}(:,3)];
    end
    min_z_val = round(min(all_z));
    max_z_val = round(max(all_z));
    
    % set up some stuff to be filled in during the plane pass
    derivs_encountered = [];
    thicks_encountered = [];
    centers_encountered = [];
    thick_means = [];
    thick_vars = [];
    deriv_means = [];
    deriv_variances = [];
    nn_dists = [];
    
    % pass the plane from min to max z value
    for l = min_z_val:max_z_val
        derivatives_here = [];
        thicks_here = [];
        centers_here = [];
        % at each plane, need to check each archaeo
        for m = 1:numel(centers_rotated_deriv)
            % how far are all points of this branch from the current plane?
            diffs_from_plane = abs(centers_rotated_deriv{m}(:,3) - l);
            [min_diff,close_ind] = min(diffs_from_plane);
            % if our closest point isn't close enough, we ignore with NaN
            if min_diff > diff_thresh
                derivatives_here(m,:) = [NaN, NaN, NaN];
                thicks_here(m,1) = NaN;
                centers_here(m,:) = [NaN, NaN, NaN];
            % if it is close enough, grab its data
            else
                derivatives_here(m,:) = outer_center_stats.derivative{m}(close_ind,:);
                thicks_here(m,1) = outer_center_stats.mean_thickness{m}(1,close_ind);
                centers_here(m,:) = outer_center_stats.spline{m}(close_ind,:);
            end
        end
        % look at this section if you want to understand change in nearest
        % neighbors
%         centers_encountered(:,:,l-min_z_val+1) = centers_here;
%         dists = squareform(pdist(centers_here));
%         dists(1:size(dists,1)+1:end) = NaN;
%         nn_dists(:,l-min_z_val+1) = min(dists,[],2);

        % set up the interesting stuff to be spit out
        thicks_encountered(:,l-min_z_val+1) = thicks_here;

        derivs_encountered(:,:,l-min_z_val+1) = derivatives_here;
        deriv_means(l-min_z_val+1,:) = mean(derivatives_here,1,'omitnan');
        deriv_variances(l-min_z_val+1,:) = var(derivatives_here,1,'omitnan');
    end 
end