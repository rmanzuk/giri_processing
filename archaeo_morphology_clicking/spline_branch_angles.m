function [branch_angles] = spline_branch_angles(branch_points_3d,center_stats,length_considered)
% This function takes the 3d branching points as well as information about
% the branch splines to calculate the angle between each branch.
%
% IN 
% branch_points_3d: n_branches x n_branches x 3 matrix where the index
% of (i,j,:) gives the 3d location of the center of the ith archaeo at the
% point of branching to the jth archaeo. 
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
    for i = 1:size(branch_points_3d,1)
        for j = 1:size(branch_points_3d,2)
            % see if these two even branch
            if branch_points_3d(i,j,1) ~= 0
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
                i_above_branch = ith_branch_spline(:,3) > branch_points_3d(i,j,3);
                i_below_thresh = i_cumlenghts < length_considered;
                i_inds = i_above_branch & i_below_thresh;
                
                j_above_branch = jth_branch_spline(:,3) > branch_points_3d(i,j,3);
                j_below_thresh = j_cumlenghts < length_considered;
                j_inds = j_above_branch & j_below_thresh;
                
                % and what is the mean heading (derivative) of each branch
                % over the interval considered
                i_mean_vec = mean(center_stats.derivative{i}(i_inds,:));
                j_mean_vec = mean(center_stats.derivative{j}(j_inds,:));
                
                % and calculate the branch angle
                branch_angles(i,j) = rad2deg(dot(i_mean_vec,j_mean_vec)/(norm(i_mean_vec)*norm(j_mean_vec)));
                
            else
            end
        end
    end
    
end
