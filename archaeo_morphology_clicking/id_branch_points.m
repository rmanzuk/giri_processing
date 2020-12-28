function [branched_flags,branching_points_3d] = id_branch_points(outers,scale_ratio,distance_threshold)
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
% branch_points_3d: n_branches x n_branches x 3 matrix where the index of
% (i,j,:) gives the 3d location of the center of the ith archaeo at the
% point of branching to the jth archaeo.

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
        for j = 1:n_branches
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
                        else
                            %do nothing because they don't branch
                        end
                    end
                end 
            end
        end
    end
end