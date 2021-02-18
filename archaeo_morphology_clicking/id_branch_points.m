function [branched_flags,branching_points_3d,outers] = id_branch_points(outers,scale_ratio,distance_threshold, combining_threshold)
% This function takes a whole slew of coral ct data and figures out which
% branches meet each other. 
%
% IN 
% outers: 1xn_branches cell array containing the cell arrays for outer circles 
% created during data collection, reshaped, but still in  
% scale_ratio: ratio of vertical image separation to pixel-width. For
% example if images are separated by 100 microns and pixels are 20 microns,
% this input is 5.
% distance_threshold: distance below which, two adjacent branches are
% supposed to come together as branches.
%
%
% OUT
% branched_flags: n_branches x n_branches comparison matrix where branches
% are kept track of. For any non-zero value in this matrix, the indices
% indicate which archaeos branch, and the values at those indices indicate
% which slice contains the branching point. Because each archaeo is
% considdered in both the rows and columns, this matrix is redundant you
% really only need the half of the data above (or below) the diagonal.
% Think of it like a covariance matrix. 
%
% branch_points_3d: n_branches x n_branches x 3 matrix where the index of
% (i,j,:) gives the 3d location of the center of the ith archaeo at the
% point of branching to the jth archaeo.
%
% outers: 1xn_branches cell array containing the cell arrays for outer circles 
% that was input to the function, but now we've accounted for branching and
% combined some that were separate before saying that they are
% continuations of the same branch.

% R. A. Manzuk 12/28/2020

    %% begin dat sweet function
    % how many slices are we working with?
    n_slices = max([cellfun(@numel,outers),cellfun(@numel,outers)],[],'all');
    
    % make sure all individual branch cell arrays are same size...just in
    % case
    for i = 1:numel(outers)
        if numel(outers{i}) < n_slices
            add_empty = cell(1,n_slices - numel(outers{i}));
            outers{i}(numel(outers{i})+1:n_slices) = add_empty;
        else
            % do nothing
        end
    end
    
    %how many branches?
    n_branches = numel(outers);
    
    % set up a flagging system to check off branches as they branch. so
    % what we're gonna do is have an n_branches x n_branches comparison
    % matrix where each branch is compared to all the others (diagonal
    % being camparison of archaeo to itself). In this matrix, the values
    % will indicate at which slice two branches come together, and the
    % indices will indicate which archaeos.
    branched_flags = zeros(n_branches);
    branching_points_3d = zeros(n_branches,n_branches,3);
    
    % start at the top of the stack and work your way down.
    for i = 1:n_slices-1
        % and then get the data for all of the unique archaeos in the slice
        slice_outers = {};
        next_outers = {};
        for j = 1:numel(outers)
            slice_outers(j) = outers{j}(i);
            next_outers(j) = outers{j}(i+1);
        end
        
        % now we need to know which ones disappear between this slice and
        % the next...maybe indicating branching
        in_this = ~cellfun(@isempty,slice_outers);
        in_next = ~cellfun(@isempty,next_outers);
        went_missing = find((in_this-in_next)==1);
        
        % will also be helpful to know which ones appear in the next slice
        will_appear = find((in_next-in_this)==1);
        
        to_trash = []; % housekeeping for branches to remove
        claimed_branches = []; % housekeeping for branches already claimed
        % if nobody (or only one) went missing, there are no branch points to worry about
        if  length(went_missing) < 2
            % do nothing
        else
            % we should go through and compare all those that went missing
            % to see if they became close and branched
            for j = went_missing
                for k = went_missing
                    if j == k
                        % do nothing
                    else
                        % calculate the distances between all points
                        distances = pdist2([slice_outers{j}(:,1),slice_outers{j}(:,2)],[slice_outers{k}(:,1),slice_outers{k}(:,2)]);
                        % if they came close enough, mark it down.
                        if min(distances(:)) < distance_threshold
                            branched_flags(j,k) = 1;
                            branched_flags(k,j) = 1;
                            branching_points_3d(j,k,1) = mean(slice_outers{j}(:,1));
                            branching_points_3d(j,k,2) = mean(slice_outers{j}(:,2));
                            branching_points_3d(j,k,3) = i * -scale_ratio;
                            branching_points_3d(k,j,1) = mean(slice_outers{k}(:,1));
                            branching_points_3d(k,j,2) = mean(slice_outers{k}(:,2));
                            branching_points_3d(k,j,3) = i * -scale_ratio;
                            
                            % and because we've recognized a branch, we
                            % should see what appears below to identify the
                            % thicker, lower branch they combine into and
                            % reallocate the points accordingly (only do
                            % this step if will_appear has values)
                            
                            % first, take the mean point of the two
                            % branches that went missing and all those that
                            % will appear
                            if ~isempty(will_appear)
                                center_missing = mean([outers{j}{i}(:,1:2);outers{k}{i}(:,1:2)]);
                                centers_appear = [];
                                count = 1;
                                for q = will_appear
                                    centers_appear(count,1:2) = mean(outers{q}{i+1}(:,1:2));
                                    count = count +1;
                                end
                                
                                %calculate the distance the center of those
                                %that went missing and each of the appearing
                                %centers
                                distances = pdist2(center_missing,centers_appear);
                                % use min of distances to index appearing
                                % branch
                                [~,min_dist_ind] = min(distances);
                                new_branch_ind = will_appear(min_dist_ind);
                                % check if this branch has been claimed
                                if ~ismember(new_branch_ind, claimed_branches)
                                    % also need to know which of the two outgoing
                                    % branches was longer (because we'll continue
                                    % that one.
                                    rel_length = [sum(~cellfun(@isempty,outers{k}(1:i))),sum(~cellfun(@isempty,outers{k}(1:i)))];
                                    [~,thicker_ind] = max(rel_length);
                                    if thicker_ind == 1
                                        to_continue = j;
                                    else
                                        to_continue = k;
                                    end
                                    
                                    % now, we just need to cut and paste and
                                    % housekeep
                                    outers{to_continue}(i+1:end) = outers{new_branch_ind}(i+1:end);
                                    to_trash = [to_trash,new_branch_ind];
                                    claimed_branches = [claimed_branches,new_branch_ind];
                                end
                            end
                        else
                            %do nothing because they don't branch
                        end
                    end
                end
            end
        end

            
        % 2 new branches appearing close togther could also indicate a branch point
        if  length(will_appear) < 2
            % do nothing
        else
            % we should go through and compare all those that will appear
            % to see if they are close and branched
            for j = will_appear
                for k = will_appear
                    if j == k
                        % do nothing
                    else
                        % calculate the distances between all points
                        distances = pdist2([next_outers{j}(:,1),next_outers{j}(:,2)],[next_outers{k}(:,1),next_outers{k}(:,2)]);
                        % if they came close enough, mark it down. But as a
                        % differnt type of brach angle (2)
                        if min(distances(:)) < distance_threshold
                            branched_flags(j,k) = 2;
                            branched_flags(k,j) = 2;
                            branching_points_3d(j,k,1) = mean(next_outers{j}(:,1));
                            branching_points_3d(j,k,2) = mean(next_outers{j}(:,2));
                            branching_points_3d(j,k,3) = (i+1) * -scale_ratio;
                            branching_points_3d(k,j,1) = mean(next_outers{k}(:,1));
                            branching_points_3d(k,j,2) = mean(next_outers{k}(:,2));
                            branching_points_3d(k,j,3) = (i+1) * -scale_ratio;
                            
                        else
                            %do nothing because they don't branch
                        end
                    end
                end
            end
        end
        % now we need to take out the trash, first looking for any
        % combinations we should have made without branching
        unique_trash = unique(to_trash);
        % which branches will appear, but haven't been combined yet
        not_assigned = setdiff(will_appear,unique_trash);
        % go through all unassigned branches and see if they should be
        % matched
        if ~isempty(went_missing)
            for k = not_assigned
                % center of the unassigned branch
                will_appear_center = mean(next_outers{k}(:,1:2));
                % centers of the gone missing branches
                went_missing_centers = [];
                count = 1;
                for j = went_missing
                    went_missing_centers(count,:) = mean(slice_outers{j}(:,1:2));
                    count = count +1;
                end
                % see if any are close enough to combine
                distances = pdist2(will_appear_center, went_missing_centers);
                [min_dist,min_dist_ind] = min(distances);
                % if close enough, combine
                if min_dist < combining_threshold
                    combine_ind = went_missing(min_dist_ind);
                    outers{combine_ind}(i+1:end) = outers{k}(i+1:end);
                    unique_trash = [unique_trash,k];
                end
            end
        end
        
        if ~isempty(to_trash)
            outers(to_trash) = [];
            outers = outers(~cellfun(@isempty,outers));

            branched_flags(to_trash,:) = [];
            branched_flags(:,to_trash) = [];

            branching_points_3d(to_trash,:,:) = [];
            branching_points_3d(:,to_trash,:) = [];
        end
        
    end
end