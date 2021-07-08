function [branch_angles, lengths_used] = spline_branch_angles2(branch_points_3d,branched_flags,center_stats, min_above_branch)
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
%
% OUT
% branch_angles: n_branches x n_branches comparison matrix where branches
% are kept track of. For any non-zero value in this matrix, the indices
% indicate which archaeos branch, and the values indicate the
% branching angle between the two.
%
% lengths_used: n_branches x n_branches comparison matrix where cell values
% are the length of branch i used to calculate the branching angle between
% branches i and j.
%
% R. A. Manzuk 02/01/2021
%% begin the function
    % loop through all possible branching combinations.
    branch_angles = zeros(size(branch_points_3d,1),size(branch_points_3d,2));
    lengths_used = zeros(size(branch_points_3d,1),size(branch_points_3d,2));
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
                
                % now, knowing the cumulative lengths and the positions
                % relative to the branching point, figure out which pieces
                % of the spline should be considdered in the calculation.
                
                % start with figuring out the two points of the branch
                % splines that are closest to each other
                spline_pdists = pdist2(ith_branch_spline,jth_branch_spline);
                [i_nearest_branch_point,j_nearest_branch_point] = find(spline_pdists==min(spline_pdists(:)));
 
                
                % account for the type of branch angle this is
                if branched_flags(i,j,1) ~= 2
                    % we know where spline i is nearest spline j (our junction)
                    % need to know ith branch's other branch points
                    i_other_branches = squeeze(branch_points_3d(i,:,:));
                    i_other_branches = i_other_branches(any(i_other_branches,2),:);
                    
                    % figure out where on ith spline those branch points
                    % are closest to.
                    i_bps_dists = pdist2(ith_branch_spline,i_other_branches);
                    [~,i_branch_spline_inds] = min(i_bps_dists);
                    
                    % and which one is next after the branch point? we'll
                    % stop our calculation there
                    if sum(i_branch_spline_inds(i_branch_spline_inds > (i_nearest_branch_point+min_above_branch))) == 0
                        i_next_ind = size(ith_branch_spline,1);
                    else
                        i_next_ind = min(i_branch_spline_inds(i_branch_spline_inds > (i_nearest_branch_point+min_above_branch)));
                    end
                    
                    % set up an array to index the proper portions of the
                    % spline
                    i_inds = false(size(ith_branch_spline,1),1);
                    i_inds(i_nearest_branch_point:i_next_ind) = true;
                    
                    % same for j
                    j_other_branches = squeeze(branch_points_3d(:,j,:));
                    j_other_branches = j_other_branches(any(j_other_branches,2),:);
        
                    j_bps_dists = pdist2(jth_branch_spline,j_other_branches);
                    [~,j_branch_spline_inds] = min(j_bps_dists);

                    if sum(j_branch_spline_inds(j_branch_spline_inds > (j_nearest_branch_point+min_above_branch))) == 0
                        j_next_ind = size(jth_branch_spline,1);
                    else
                        j_next_ind = min(j_branch_spline_inds(j_branch_spline_inds > (j_nearest_branch_point+min_above_branch)));
                    end
                   
                    j_inds = false(size(jth_branch_spline,1),1);
                    j_inds(j_nearest_branch_point:j_next_ind) = true;
                    
                elseif branched_flags(i,j,1) == 2
                    % for this kind, consider enties before and up to branch point 
                    i_other_branches = squeeze(branch_points_3d(i,:,:));
                    i_other_branches = i_other_branches(any(i_other_branches,2),:);
                    
                    % figure out where on ith spline those branch points
                    % are closest to.
                    i_bps_dists = pdist2(ith_branch_spline,i_other_branches);
                    [~,i_branch_spline_inds] = min(i_bps_dists);
                    
                    % and which one is next after the branch point? we'll
                    % stop our calculation there
                    if sum(i_branch_spline_inds(i_branch_spline_inds < (i_nearest_branch_point-min_above_branch))) == 0
                        i_next_ind = 1;
                    else
                        i_next_ind = max(i_branch_spline_inds(i_branch_spline_inds < (i_nearest_branch_point-min_above_branch)));
                    end
                    
                    % set up an array to index the proper portions of the
                    % spline
                    i_inds = false(size(ith_branch_spline,1),1);
                    i_inds(i_next_ind:i_nearest_branch_point) = true;
                    
                    % same for j
                    j_other_branches = squeeze(branch_points_3d(:,j,:));
                    j_other_branches = j_other_branches(any(j_other_branches,2),:);
        
                    j_bps_dists = pdist2(jth_branch_spline,j_other_branches);
                    [~,j_branch_spline_inds] = min(j_bps_dists);

                    if sum(j_branch_spline_inds(j_branch_spline_inds < (j_nearest_branch_point-min_above_branch))) == 0
                        j_next_ind = 1;
                    else
                        j_next_ind = max(j_branch_spline_inds(j_branch_spline_inds < (j_nearest_branch_point-min_above_branch)));
                    end
                   
                    j_inds = false(size(jth_branch_spline,1),1);
                    j_inds(j_next_ind:j_nearest_branch_point) = true;
                    
                end
                % lets make vectors between all points and the first point
                % of the parts of the spline considered
                if branched_flags(i,j,1) ~= 2
                    i_considered = ith_branch_spline(i_inds,:);
                    i_all_vecs = i_considered(2:end,:) - i_considered(1,:);
                    j_considered = jth_branch_spline(j_inds,:);
                    j_all_vecs = j_considered(2:end,:) - j_considered(1,:);
                elseif branched_flags(i,j,1) == 2
                    i_considered = ith_branch_spline(i_inds,:);
                    i_all_vecs = i_considered(1:(end-1),:) - i_considered(end,:);
                    j_considered = jth_branch_spline(j_inds,:);
                    j_all_vecs = j_considered(1:(end-1),:) - j_considered(end,:);
                end
                
                % make those unit vectors and take their average
                i_norms = vecnorm(i_all_vecs')';
                j_norms = vecnorm(j_all_vecs')';
                i_unit_vecs = i_all_vecs./i_norms;
                j_unit_vecs = j_all_vecs./j_norms;
                i_vec = mean(i_unit_vecs);
                j_vec = mean(j_unit_vecs);

                % and calculate the branch angle
                if numel(i_vec) == 3 && numel(j_vec) == 3
                    branch_angles(i,j) = rad2deg(atan2(norm(cross(i_vec,j_vec)),dot(i_vec,j_vec)));
                else
                    branch_angles(i,j) = 0;
                end
                
                % lastly, we can calculate the length of each branch used in the calculation. 
                if length(center_stats.spline{i}(i_inds,1)) > 2
                    x = center_stats.spline{i}(i_inds,1);
                    y = center_stats.spline{i}(i_inds,2);
                    z = center_stats.spline{i}(i_inds,3);
                
                    lengths_used(i,j) = arclength(x,y,z);
                else
                    lengths_used(i,j) = 0;
                end
                    
            else
            end
        end
    end
    
end
