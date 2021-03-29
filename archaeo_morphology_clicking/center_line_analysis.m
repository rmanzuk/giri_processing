function[center_stats] = center_line_analysis(outlines_3d,center_lines,sample_freq,thickness_sampling)
% This function takes 3d archaeo point clouds (densified or not) finds the
% center line of each archaeo at the disired sampling resolution, with the
% point clouds as is with no rotations.
%
% IN
% outlines_3d: 1xn_branches cell array containing the densified or non-densified 3d outputs for
% inner clicked data from the densify_3d or make_clicking_3d functions.
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
% center_stats: Structure containing 1xn_branch cell arrays where each 
% cell contains a given measure running along that branch.
% These measurements in order are:[in
%   - spline:the point values from evaluation of a spline fit to the center data
%   - derivative: the point values from evaluation of derivative of the spline.
%   - inclinations: local inclination at every evaluation point.
%   - min_thickness: distance from centroid to nearest point at each evaluation interval.
%   - max_thickness: distance from centroid to furthest point at each evaluation interval.
%   - mean_thicknes: mean distance from centroid to all points at each evaluation interval.
%   - new_centroids: positions of centers resulting from thickness analysis.
%
%
% Ryan A. Manzuk 09/08/2020
    %% begin the function
    center_spline_points = {};
    derivative_spline_points = {};
    inclinations = {};
    declinations = {};
    min_thicknesses = {};
    max_thicknesses = {};
    mean_thicknesses = {};
    new_centroids = {};
    
    % first we actually need to know which center lines are long enough to
    % analyze
    long_enough = logical(zeros(1,numel(center_lines)));
    
    for k = 1:numel(center_lines)
        xyz = center_lines{k};
        xyz = xyz(all(~isnan(xyz),2),:);
        if numel(xyz) > 6
            long_enough(k) = true;
        end
    end

    for i = find(long_enough)
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
        [declines,slope_run] = cart2pol(derivative_eval(:,1),derivative_eval(:,2));
        inclines = atand(derivative_eval(:,3)./slope_run);
        declines = rad2deg(declines);
        declines = -1.*(declines - 90);
        declines(declines <=0) = 360 + declines(declines <=0);

        % now based upon local curves, 
        min_thick = [];
        max_thick = [];
        mean_thick = [];
        center_points = [];
        for j = 1:n_eval
            if ~isempty(outlines_3d{i})
                % make the centerline spline pointthe middle of the point cloud. 
                new_archaeo_cloud = outlines_3d{i}(:,1:3) - center_line_eval(j,:);
                % rotate the cloud based on the local derivative
                % first, rotation matrix for spin around z (using declination)axis such that the
                % tangent is facing north.
                z_rot_mat = [cosd(declines(j)), -sind(declines(j)), 0;
                sind(declines(j)), cosd(declines(j)), 0; 0, 0, 1];
                % then rotation matrix to account for inclinations
                x_rot_mat = [1,0,0;0,cosd(90-inclines(j)), -sind(90-inclines(j));
                0, sind(90-inclines(j)), cosd(90-inclines(j))];
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
            else
                min_thick(j) = NaN;
                max_thick(j) = NaN;
                mean_thick(j) = NaN; 
                center_points(j,:) = [NaN, NaN, NaN];
            end
        end
        % set all of the statistics for this archaeo
        center_spline_points{i} = center_line_eval;
        derivative_spline_points{i} = derivative_eval;
        inclinations{i} = inclines;
        declinations{i} = declines;
        min_thicknesses{i} = min_thick;
        max_thicknesses{i} = max_thick;
        mean_thicknesses{i} = mean_thick;
        new_centroids{i} = center_points;
    end
    
    % and just need to ensure sizes are kept, so put empty arrays where not
    % long enough
    for i = find(~long_enough)
        center_spline_points{i} = [];
        derivative_spline_points{i} = [];
        inclinations{i} = [];
        declinations{i} = [];
        min_thicknesses{i} = [];
        max_thicknesses{i} = [];
        mean_thicknesses{i} = [];
        new_centroids{i} = [];
    end
    
    % and set up structures for easy output
    center_stats.spline = center_spline_points;
    center_stats.derivative = derivative_spline_points;
    center_stats.inclinations = inclinations;
    center_stats.declinations = declinations;
    center_stats.min_thickness = min_thicknesses;
    center_stats.max_thickness = max_thicknesses;
    center_stats.mean_thickness = mean_thicknesses;
    center_stats.new_centroids = new_centroids;
    
end
