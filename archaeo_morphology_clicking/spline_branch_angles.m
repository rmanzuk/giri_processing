function [branch_angles, radii_curvature] = spline_branch_angles(branch_points_3d,branched_flags,center_stats,length_considered)
% This function takes the 3d branching points as well as information about
% the branch splines to calculate the angle between each branch.
%
% IN 
% branch_points_3d: n_branches x n_branches x 3 matrix where the index
% of (i,j,:) gives the 3d location of the center of the ith archaeo at the
% point of branching to the jth archaeo. 
%
% branched_flags: n_branches x n_branches comparison matrix where branches
% are kept track of. For any non-zero value in this matrix, the indices
% indicate which archaeos branch, and the values at those indices indicate
% which slice contains the branching point. Because each archaeo is
% considdered in both the rows and columns, this matrix is redundant you
% really only need the half of the data above (or below) the diagonal.
% Think of it like a covariance matrix. 
%
% center_stats: Structure containing
% 1xn_branch cell arrays where each cell contains a given measure of the
% branch. Output from center_line_analysis.m 
%
% length_considered: interval
% (in the units of the coordinate space) above the branch point that should
% be considered in the branch angle measurement. So if your coordinate
% space is in microns, to measure branches based upon 100 microns above the
% branch point, this input should be 100.
%
%
% OUT
% branch_angles: n_branches x n_branches comparison matrix where branches
% are kept track of. For any non-zero value in this matrix, the indices
% indicate which archaeos branch, and the values indicate the
% branching angle between the two.
%
% R. A. Manzuk 02/01/2021
%% begin the function
    % loop through all possible branching combinations.
    branch_angles = zeros(size(branch_points_3d,1),size(branch_points_3d,2));
    radii_curvature = zeros(size(branch_points_3d,1),size(branch_points_3d,2));
    for i = 1:size(branch_points_3d,1)
        for j = 1:size(branch_points_3d,2)
            % see if these two even branch and are long enough
            if branch_points_3d(i,j,1) ~= 0 && ~isempty(center_stats.spline{i}) ...
                    && ~isempty(center_stats.spline{j}) && ...
                    sum(~isnan(center_stats.spline{i}(:,1))) > 1 && ...
                    sum(~isnan(center_stats.spline{j}(:,1))) > 1 
                
                % extract the splines for the two branches
                ith_branch_spline = center_stats.spline{i};
                jth_branch_spline = center_stats.spline{j};
                
                % and figure out the progressive length of the spline
                [~,i_lengths] = arclength(ith_branch_spline(:,1),ith_branch_spline(:,2),ith_branch_spline(:,3));
                i_cumlenghts = [0;cumsum(i_lengths)];
                [~,j_lengths] = arclength(jth_branch_spline(:,1),jth_branch_spline(:,2),jth_branch_spline(:,3));
                j_cumlenghts = [0;cumsum(j_lengths)];
                
                % now, knowing the cumulative lengths and the positions
                % relative to the branching point, figure out which pieces
                % of the spline should be considdered in the calculation.
                [~,i_nearest_branch_point] = min(abs(ith_branch_spline(:,3) - branch_points_3d(i,j,3)));
                i_past_branch = false(size(ith_branch_spline,1),1);
                
                [~,j_nearest_branch_point] = min(abs(jth_branch_spline(:,3) - branch_points_3d(i,j,3)));
                j_past_branch = false(size(jth_branch_spline,1),1);
                
                % account for the type of branch angle this is
                if branched_flags(i,j,1) ~= 2
                    % for this kind, can consider entries after branch
                    % point 
                    i_past_branch(i_nearest_branch_point:end) = true;
                    % then need cum_lengths from branch point
                    i_new_cumlengths = i_cumlenghts - i_cumlenghts(i_nearest_branch_point,:);
                    i_below_thresh = i_new_cumlengths < length_considered;
                    i_inds = i_past_branch & i_below_thresh;
                    
                    % same for j
                    j_past_branch(j_nearest_branch_point:end) = true;
                    % then need cum_lengths from branch point
                    j_new_cumlengths = j_cumlenghts - j_cumlenghts(j_nearest_branch_point,:);
                    j_below_thresh = j_new_cumlengths < length_considered;
                    j_inds = j_past_branch & j_below_thresh;
                    
                elseif branched_flags(i,j,1) == 2
                    % for this kind, consider enties before and up to branch point 
                    i_past_branch(1:i_nearest_branch_point) = true;
                    % then need cum_lengths from branch point
                    i_new_cumlengths = i_cumlenghts - i_cumlenghts(i_nearest_branch_point,:);
                    i_below_thresh = i_new_cumlengths > -length_considered;
                    i_inds = i_past_branch & i_below_thresh;
                    
                    % same for j
                    j_past_branch(1:j_nearest_branch_point) = true;
                    % then need cum_lengths from branch point
                    j_new_cumlengths = j_cumlenghts - j_cumlenghts(j_nearest_branch_point,:);
                    j_below_thresh = j_new_cumlengths >- length_considered;
                    j_inds = j_past_branch & j_below_thresh;
                end
                
                
                % and what is the mean heading (derivative) of each branch
                % over the interval considered
                i_mean_vec = mean(center_stats.derivative{i}(i_inds,:));
                j_mean_vec = mean(center_stats.derivative{j}(j_inds,:));

                % and calculate the branch angle
                if isequal(size(i_mean_vec), size(j_mean_vec))
                    branch_angles(i,j) = abs(acosd(dot(i_mean_vec,j_mean_vec)/(norm(i_mean_vec)*norm(j_mean_vec))));
                else
                    branch_angles(i,j) = 0;
                end
                
                % lastly, we can calculate the radius of curvature for the
                % ith branch at this branch point. 
                x = center_stats.spline{i}(i_inds,1);
                y = center_stats.spline{i}(i_inds,2);
                z = center_stats.spline{i}(i_inds,3);
                
                triangulation = delaunayTriangulation(x,y,z);
                [~,r] = circumcenter(triangulation);
                radii_curvature(i,j) = mean(r);
            else
            end
        end
    end
    
end
