function[inner_center_stats,outer_center_stats] = center_line_analysis(inner_3d,outer_3d,center_lines,sample_freq,thickness_sampling)
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
% center_lines: 1xn_archaeo cell array where each cell contains the 3D
% coordinates for the center line of the outer tube for the given archaeo.
%
% sampling_freq: interval (in the units of the coordinate space) at
% which you would like to sample the center lines for statistics. So if your
% coordinate space is in microns, and you want to sample every 5 microns,
% this variable should be 5.
%
% thickness_sampling: interval (in the units of the coordinate space) above
% and below a center point that should be included when evaluating
% thickness. If this variable is set to 1, then 1 unit above and below the
% point will be considdered in thickness measurements.
% 
% OUT
% inner_center_stats: Structure containing 1xn_archaeo cell arrays where each 
% cell contains a given measure of the inner tube for the given archaeo.
% These measurements in order are:[in
%   - spline:the point values from evaluation of a spline fit to the center data
%   - derivative: the point values from evaluation of derivative of the spline.
%   - inclinations: local inclination at every evaluation point.
%   - min_thickness: distance from centroid to nearest point at each evaluation interval.
%   - max_thickness: distance from centroid to furthest point at each evaluation interval.
%   - mean_thicknes: mean distance from centroid to all points at each evaluation interval.
%   - new_centroids: positions of centers resulting from thickness analysis.
%
% outer_center_stats: Structure containing 1xn_archaeo cell arrays where each 
% cell contains a given measure of the outer tube for the given archaeo.
% These measurements in order are:
%   - spline:the point values from evaluation of a spline fit to the center data
%   - derivative: the point values from evaluation of derivative of the spline.
%   - inclinations: local inclination at every evaluation point.
%   - min_thickness: distance from centroid to nearest point at each evaluation interval.
%   - max_thickness: distance from centroid to furthest point at each evaluation interval.
%   - mean_thicknes: mean distance from centroid to all points at each evaluation interval.
%   - new_centroids: positions of centers resulting from thickness analysis.
%
% Ryan A. Manzuk 09/08/2020
    %% begin the function
    inner_center_spline_points = {};
    inner_derivative_spline_points = {};
    inner_inclinations = {};
    inner_declinations = {};
    inner_min_thicknesses = {};
    inner_max_thicknesses = {};
    inner_mean_thicknesses = {};
    inner_new_centroids = {};

    for i = 1:numel(center_lines)
        % extract individual archaeo center line and get rid of nans
        xyz = center_lines{i}';
        xyz = xyz(:,all(~isnan(xyz)));
        % assess center line arclength so we know how many samples to take of
        % its spline
        curve_length = arclength(xyz(1,:),xyz(2,:),xyz(3,:));
        n_eval = round(curve_length/sample_freq);
        % spline fit of the center line and its derivative
        center_spline = cscvn(xyz);
        spline_derivative = fnder(center_spline);
        % set up a set of values to evaluate and do so on the spline and its
        % derivative
        to_eval = linspace(0,center_spline.breaks(end), n_eval);
        center_line_eval = ppval(center_spline,to_eval)';
        derivative_eval = ppval(spline_derivative,to_eval)';

        %turn derivative at each point into inclination and declination
        [declinations,slope_run] = cart2pol(derivative_eval(:,1),derivative_eval(:,2));
        inclinations = atand(derivative_eval(:,3)./slope_run);
        declinations = rad2deg(declinations);
        declinations = -1.*(declinations - 90);
        declinations(declinations <=0) = 360 + declinations(declinations <=0);

        % now based upon local curves, 
        min_thick = [];
        max_thick = [];
        mean_thick = [];
        center_points = [];
        for j = 1:n_eval
            % make the centerline spline pointthe middle of the point cloud. 
            new_archaeo_cloud = inner_3d{i}(:,1:3) - center_line_eval(j,:);
            % rotate the cloud based on the local derivative
            % first, rotation matrix for spin around z (using declination)axis such that the
            % tangent is facing north.
            z_rot_mat = [cosd(declinations(j)), -sind(declinations(j)), 0;
            sind(declinations(j)), cosd(declinations(j)), 0; 0, 0, 1];
            % then rotation matrix to account for inclinations
            x_rot_mat = [1,0,0;0,cosd(90-inclinations(j)), -sind(90-inclinations(j));
            0, sind(90-inclinations(j)), cosd(90-inclinations(j))];
            % rotate it!
            rotated_1 = z_rot_mat * new_archaeo_cloud';
            rotated_2 = (x_rot_mat * rotated_1)';
            % take the points above and below the proper z values
            points_in_window = rotated_2(:,3) < thickness_sampling & rotated_2(:,3) > -thickness_sampling;
            if sum(points_in_window) >= 10
                this_chunk = rotated_2(points_in_window,:);
                % what's the thickness?
                centroid = mean(this_chunk(:,1:2));
                distances = sqrt((centroid(1)-this_chunk(:,1)).^2 + (centroid(2)-this_chunk(:,2)).^2);
                min_thick(j) = min(distances, [],'all');
                max_thick(j) = max(distances, [],'all');
                mean_thick(j) = mean(distances,'all');
                center_points(j,:) = (-z_rot_mat*(-x_rot_mat*[centroid,0]'))' + center_line_eval(j,:);
            else
                % if you don't have enough points, don't try to do anything
                min_thick(j) = NaN;
                max_thick(j) = NaN;
                mean_thick(j) = NaN; 
                center_points(j,:) = [NaN, NaN, NaN];
            end
        end
        % set all of the statistics for this archaeo
        inner_center_spline_points{i} = center_line_eval;
        inner_derivative_spline_points{i} = derivative_eval;
        inner_inclinations{i} = inclinations;
        inner_declinations{i} = declinations;
        inner_min_thicknesses{i} = min_thick;
        inner_max_thicknesses{i} = max_thick;
        inner_mean_thicknesses{i} = mean_thick;
        inner_new_centroids{i} = center_points;
    end

    outer_center_spline_points = {};
    outer_derivative_spline_points = {};
    outer_inclinations = {};
    outer_min_thicknesses = {};
    outer_max_thicknesses = {};
    outer_mean_thicknesses = {};
    outer_new_centroids = {};

    for i = 1:numel(center_lines)
        % extract individual archaeo center line and get rid of nans
        xyz = center_lines{i}';
        xyz = xyz(:,all(~isnan(xyz)));
        % assess center line arclength so we know how many samples to take of
        % its spline
        curve_length = arclength(xyz(1,:),xyz(2,:),xyz(3,:));
        n_eval = round(curve_length/sample_freq);
        % spline fit of the center line and its derivative
        center_spline = cscvn(xyz);
        spline_derivative = fnder(center_spline);
        % set up a set of values to evaluate and do so on the spline and its
        % derivative
        to_eval = linspace(0,center_spline.breaks(end), n_eval);
        center_line_eval = ppval(center_spline,to_eval)';
        derivative_eval = ppval(spline_derivative,to_eval)';

        %turn derivative at each point into inclination and declination
        [declinations,slope_run] = cart2pol(derivative_eval(:,1),derivative_eval(:,2));
        inclinations = atand(derivative_eval(:,3)./slope_run);
        declinations = rad2deg(declinations);
        declinations = -1.*(declinations - 90);
        declinations(declinations <=0) = 360 + declinations(declinations <=0);

        % now based upon local curves, 
        min_thick = [];
        max_thick = [];
        mean_thick = [];
        center_points = [];
        for j = 1:n_eval
            % make the centerline spline pointthe middle of the point cloud. 
            new_archaeo_cloud = outer_3d{i}(:,1:3) - center_line_eval(j,:);
            % rotate the cloud based on the local derivative
            % first, rotation matrix for spin around z axis such that the
            % derivative is facing north.
            z_rot_mat = [cosd(declinations(j)), -sind(declinations(j)), 0;
            sind(declinations(j)), cosd(declinations(j)), 0; 0, 0, 1];
            % then rotation matrix to account for inclinations
            x_rot_mat = [1,0,0;0,cosd(90-inclinations(j)), -sind(90-inclinations(j));
            0, sind(90-inclinations(j)), cosd(90-inclinations(j))];
            % rotate it!
            rotated_1 = z_rot_mat * new_archaeo_cloud';
            rotated_2 = (x_rot_mat * rotated_1)';
            % take the points above and below the proper z values
            points_in_window = rotated_2(:,3) < thickness_sampling & rotated_2(:,3) > -thickness_sampling;
            if sum(points_in_window) >= 10
                this_chunk = rotated_2(points_in_window,:);
                % what's the thickness?
                centroid = mean(this_chunk(:,1:2));
                distances = sqrt((centroid(1)-this_chunk(:,1)).^2 + (centroid(2)-this_chunk(:,2)).^2);
                min_thick(j) = min(distances, [],'all');
                max_thick(j) = max(distances, [],'all');
                mean_thick(j) = mean(distances,'all');
                center_points(j,:) = (z_rot_mat'*(x_rot_mat'*[centroid,0]'))' + center_line_eval(j,:);
            else
                min_thick(j) = NaN;
                max_thick(j) = NaN;
                mean_thick(j) = NaN; 
                center_points(j,:) = [NaN, NaN, NaN];
            end
        end
        outer_center_spline_points{i} = center_line_eval;
        outer_derivative_spline_points{i} = derivative_eval;
        outer_inclinations{i} = inclinations;
        outer_declinations{i} = declinations;
        outer_min_thicknesses{i} = min_thick;
        outer_max_thicknesses{i} = max_thick;
        outer_mean_thicknesses{i} = mean_thick;
        outer_new_centroids{i} = center_points;
    end
    
    % and set up structures for easy output
    inner_center_stats.spline = inner_center_spline_points;
    inner_center_stats.derivative = inner_derivative_spline_points;
    inner_center_stats.inclinations = inner_inclinations;
    inner_center_stats.declinations = inner_declinations;
    inner_center_stats.min_thickness = inner_min_thicknesses;
    inner_center_stats.max_thickness = inner_max_thicknesses;
    inner_center_stats.mean_thickness = inner_mean_thicknesses;
    inner_center_stats.new_centroids = inner_new_centroids;
    
    outer_center_stats.spline = outer_center_spline_points;
    outer_center_stats.derivative = outer_derivative_spline_points;
    outer_center_stats.inclinations = outer_inclinations;
    outer_center_stats.declinations = outer_declinations;
    outer_center_stats.min_thickness = outer_min_thicknesses;
    outer_center_stats.max_thickness = outer_max_thicknesses;
    outer_center_stats.mean_thickness = outer_mean_thicknesses;
    outer_center_stats.new_centroids = outer_new_centroids;
end
