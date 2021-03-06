function[center_points] = iterate_center_lines(outer_3d,points_here_thresh,iterate_stop,n_iter,sample_freq,thickness_sampling,derivative_smooth)
% This function takes 3d archaeo point clouds (densified or not) finds the
% center line of each archaeo at the disired sampling resolution, with the
% point clouds as is with no rotations. It does so iteratively, moving the
% center line each time based upon a new centroid informed by the local
% heading of the archeo.
%
% IN
%
% outer_3d: 1xn_archaeos cell array containing the densified or non-densified 3d outputs for
% outer clicked data from the densify_3d or make_clicking_3d functions.
%
% points_here_thresh: minimum number of points that need to be within an
% interval for the function to assess the centroid for the interval.
%
% iterate_stop: Threshold for similarity between iterations. If The
% difference between the centerline matrix between 2 iterations is below
% this threshold, the loop will stop iterating and give the latest
% iteration as the final center line.
%
% n_iter: number of iterations to go through if iterate_stop never is
% activated
%
% sample_freq: interval (in the units of the coordinate space) at
% which you would like to sample the center lines for statistics. So if your
% coordinate space is in microns, and you want to sample every 5 microns,
% this variable should be 5.
%
% thickness_sampling: interval (in the units of the coordinate space) above
% and below a center point that should be included when evaluating
% thickness. If this variable is set to 1, then 1 unit above and below the
% point will be considdered in thickness measurements.
%
% derivative_smooth: number of points for the moving mean of the derivative
% to help not take wild twists into account
% 
% OUT
% center_points: 1xn_branch cell array where each cell contains the 3D
% coordinates for the center line of the given branch.
%
% Ryan A. Manzuk 09/22/2020

    %% begin the function
    center_points = {};
    for i = 1:numel(outer_3d)
        % extract 3d (dense) data for this archaeo
        current_points = outer_3d{i};
        % set up xyz data matrix
        xyz = current_points(:,1:3);

        % get a starting center line based upon the branch being rotated to
        % its best orientation.
        init_centers = svd_center_lines({xyz},thickness_sampling,points_here_thresh);
        % now just set xyz to be the center points and remove nans.
        xyz = init_centers{1}';
        xyz = xyz(:,all(~isnan(xyz)));
        % we'll now need to evaluate a spline going through our initial
        % center points. Set that up for a number of evaluations given
        % the length of the curve through the inial centers
        if size(xyz,2) > 2
            curve_length = arclength(xyz(1,:),xyz(2,:),xyz(3,:));
            n_eval = round(curve_length/sample_freq);
            for k = 1:n_iter
                % spline fit of the center line and its derivative
                center_spline = cscvn(xyz);
                spline_derivative = fnder(center_spline);
                % set up a set of values to evaluate and do so on the spline and its
                % derivative
                to_eval = linspace(0,center_spline.breaks(end), n_eval);
                center_line_eval = ppval(center_spline,to_eval)';
                derivative_eval = ppval(spline_derivative,to_eval)';
                
                % from the derivative, take the mean, so we know the mean
                % direction of movement and can decide the best axes to
                % smooth
                mean_direction = mean(derivative_eval);
                x_axis_angle = acosd(dot([1,0,0],mean_direction)/(norm([1,0,0])*norm(mean_direction)));
                y_axis_angle = acosd(dot([0,1,0],mean_direction)/(norm([0,1,0])*norm(mean_direction)));
                z_axis_angle = acosd(dot([0,0,1],mean_direction)/(norm([0,0,1])*norm(mean_direction)));
                
                % which angle is smallest? Don't smooth that 
                init_centers = [];
                if x_axis_angle < y_axis_angle && x_axis_angle < z_axis_angle
                    init_centers = [center_line_eval(:,1),movmean(center_line_eval(:,2),derivative_smooth),movmean(center_line_eval(:,3),derivative_smooth)];
                elseif y_axis_angle < x_axis_angle && y_axis_angle < z_axis_angle
                    init_centers = [movmean(center_line_eval(:,1),derivative_smooth), center_line_eval(:,2), movmean(center_line_eval(:,3),derivative_smooth)];
                else
                    init_centers = [movmean(center_line_eval(:,1),derivative_smooth), movmean(center_line_eval(:,2),derivative_smooth), center_line_eval(:,3)];
                end
                
                % then redo the center spline and derivative with the
                % smoother line
                center_spline = cscvn(init_centers');
                spline_derivative = fnder(center_spline);
                to_eval = linspace(0,center_spline.breaks(end), n_eval);
                center_line_eval = ppval(center_spline,to_eval)';
                derivative_eval = ppval(spline_derivative,to_eval)';
                
                % and if any parts of the derivative eval are outliers,
                % remove them and replace with the mean direction
                [~,outlying_derivatives] = rmoutliers(derivative_eval);
                derivative_eval(outlying_derivatives,:) = repmat(mean_direction,sum(outlying_derivatives),1);
                
                
                %turn derivative at each point into inclination and declination
                [declinations,slope_run] = cart2pol(derivative_eval(:,1),derivative_eval(:,2));
                inclinations = atand(derivative_eval(:,3)./slope_run);
                declinations = rad2deg(declinations);
                declinations = -1.*(declinations - 90);
                declinations(declinations <=0) = 360 + declinations(declinations <=0);
                
                % to smooth the derivative and account for fewer wild
                % turns, first we need to know which axis to keep
                % unsmoothed 

                new_center_points = [];
                for l = 1:size(center_line_eval,1)
                    % make the centerline spline pointthe middle of the point cloud. 
                    new_archaeo_cloud = outer_3d{i}(:,1:3) - center_line_eval(l,:);
                    % rotate the cloud based on the local derivative
                    % first, rotation matrix for spin around z (using declination)axis such that the
                    % tangent is facing north.
                    z_rot_mat = [cosd(declinations(l)), -sind(declinations(l)), 0;
                    sind(declinations(l)), cosd(declinations(l)), 0; 0, 0, 1];
                    % then rotation matrix to account for inclinations
                    x_rot_mat = [1,0,0;0,cosd(90-inclinations(l)), -sind(90-inclinations(l));
                    0, sind(90-inclinations(l)), cosd(90-inclinations(l))];
                    % rotate it!
                    rotated_1 = z_rot_mat * new_archaeo_cloud';
                    rotated_2 = (x_rot_mat * rotated_1)';
                    % take the points above and below the proper z values
                    points_in_window = rotated_2(:,3) < thickness_sampling & rotated_2(:,3) > -thickness_sampling;
                    if sum(points_in_window) >= 10
                        this_chunk = rotated_2(points_in_window,:);
                        % where is the new centroid
                        centroid = mean(this_chunk);
                        new_center_points(l,:) = (z_rot_mat'*(x_rot_mat'*centroid'))' + center_line_eval(l,:);
                    else
                        % if you don't have enough points, don't try to do anything
                        new_center_points(l,:) = [NaN, NaN, NaN];
                    end
                end
                % now we need to compare the new center line and the old
                % center line...if the new one is not that different from
                % the old one... break the loop
                difference_mat = new_center_points(:,:) - center_line_eval;
                squared_diff = difference_mat.^2;
                not_nan = ~isnan(squared_diff);
                final_diff_metric = sum(squared_diff(not_nan),'all');
                if final_diff_metric < iterate_stop
                    new_center_points = rmoutliers(new_center_points);
                    center_points{i} = new_center_points(:,:);
                    break
                else
                    if k == n_iter
                        new_center_points = rmoutliers(new_center_points);
                        center_points{i} = new_center_points(:,:);
                        break
                    end
                    xyz = new_center_points(:,:)';
                    xyz = xyz(:,all(~isnan(xyz)));
                end
            end 
        else
            center_points{i} = [NaN,NaN,NaN];
        end
    end
end