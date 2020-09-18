function [branched_flags,branching_angles,branching_points_3d] = process_branched_network(inners,outers,scale_ratio,slices_above)
% This function takes a whole slew of archaeo data and figures out which
% archaeos branch to each other and what the angle of branching is. Right
% now, the the indicator of branching is having overlapping outer circle
% datapoints. That method could easily be changed if another were thought
% up. 
%
% IN 
% inners: 1xn_archaeos cell array containing the cell arrays for inner circles 
% created during data collection. To set this variable up, I just call 
% inners = {inner1, inner2, ..., innern}
% outers: 1xn_archaeos cell array containing the cell arrays for outer circles 
% created during data collection. To set this variable up, I just call 
% outers = {outer1, outer2, ..., outern}
% scale_ratio: ratio of vertical image separation to pixel-width. For
% example if images are separated by 100 microns and pixels are 20 microns,
% this input is 5.
% slices_above: the number of slices above a branching point which will be
% taken into account when calculating the direction vectors for branching
% archaeos.
%
% OUT
% branched_flags: n_archaeos x n_archaeos comparison matrix where branches
% are kept track of. For any non-zero value in this matrix, the indices
% indicate which archaeos branch, and the values at those indices indicate
% which slice contains the branching point. Because each archaeo is
% considdered in both the rows and columns, this matrix is redundant you
% really only need the half of the data above (or below) the diagonal.
% Think of it like a covariance matrix. 
% branching_angles: same as branched_flags but the values indicate the
% branching angle between the two archaeos as opposed to the branch point
% slice number.
% branch_points_3d: n_archaeos x n_archaeos x 3 matrix where the index of
% (i,j,:) gives the 3d location of the center of the ith archaeo at the
% point of branching to the jth archaeo.

% R. A. Manzuk 07/18/2020

    %% begin dat sweet function
    % how many slices are we working with?
    n_slices = max([cellfun(@numel,inners),cellfun(@numel,outers)],[],'all');

    % make sure all individual archaeo cell arrays are same size...just in
    % case
    for i = 1:numel(inners)
        if numel(inners{i}) < n_slices
            add_empty = cell(1,n_slices - numel(inners{i}));
            inners{i}(numel(inners{i})+1:n_slices) = add_empty;
        else
            % do nothing
        end
    end

    for i = 1:numel(outers)
        if numel(outers{i}) < n_slices
            add_empty = cell(1,n_slices - numel(outers{i}));
            outers{i}(numel(outers{i})+1:n_slices) = add_empty;
        else
            % do nothing
        end
    end

    %how many archaeos?
    n_archaeos = numel(inners);

    % set up a flagging system to check off archaeos as they branch. so
    % what we're gonna do is have an n_archaeos x n_archaeos comparison
    % matrix where each archaeo is compared to all the others (diagonal
    % being camparison of archaeo to itself). In this matrix, the values
    % will indicate at which slice two archaeos come together, and the
    % indices will indicate which archaeos.
    branched_flags = zeros(n_archaeos);
    
    % start at the top of the stack and work your way down.
    for i = 1:n_slices    
        % and then get the data for all of the unique archaeos in the slice
        slice_inners = {};
        slice_outers = {};
        for j = 1:n_archaeos
            slice_inners(j) = inners{j}(i);
            slice_outers(j) = outers{j}(i);
        end

        % now we need to know which ones overlap aka branch
        warning('off','all');
        pgons = [];
        for j = 1:n_archaeos
            % define the polygon for each archaeo in the slice using outers
            % data
            if isempty(slice_outers{j})
                poly = polyshape();
            else
                poly = polyshape(slice_outers{j}(:,1),slice_outers{j}(:,2));
            end
            pgons = [pgons, poly];
        end
        warning('on','all');
           
        % if two archaeo polygons overlap, we'll call that a branch
        touchers = overlaps(pgons);
        
        if sum(touchers) == 0
            % do nothing because there are no branches in this slice
        else
            % get rid of the diagonal values becuase of course each archaeo
            % is going to overlap with itself
            touchers = touchers .*  ~diag(ones(1,n_archaeos));
            % if this is the first time in the stack where we've recognized
            % any branches, we don't need to worry about redundancy, so we
            % can just assign these branches to our branched flags (multiplied by i to indicate
            % slice number) 
            if sum(branched_flags) == 0
                branched_flags = touchers.*i;
            % otherwise, we need to check for redundancy so we don't
            % overwrite branches. In this system, we would never double
            % count a branch....but we want to maintain the first slice
            % that a branch is recognized as the datapoint and ignore that
            % branch later down
            else
                % we account for redundancy by just subtracting the indices
                % of the branches already recognized
                new_branches = (touchers - (branched_flags~=0));
                % and multiply by i to indicate slice number
                new_branches2 = (new_branches==1).*i;
                % and add the new branches to our existing flagging matrix
                branched_flags = branched_flags + new_branches2;
            end
        end
    end

    % cool, we know which archaeos branch to each other, and in which slice
    % number this happens. Now, we just calculate the angles.
    
    % Keep track of angles in a similar comparison matrix where in stead of
    % noting slice number of the branching point, we'll note angle.
    branching_angles = zeros(n_archaeos);
    branching_points_3d = zeros(n_archaeos,n_archaeos,3);
    
    % Go through each entry of the flagging matrix
    for i = 1:n_archaeos
        for j = 1:n_archaeos
            % if it's not zero, we know we've got a branch
            if branched_flags(i,j) ~= 0
                % and take the ranch of slices from which we'll calculate
                % the vectors of each archaeo for branch angle
                branch_slice = branched_flags(i,j);
                top_slice = branch_slice-slices_above+1;
                % can't have a top slice outside the stack
                if top_slice <1
                    top_slice = 1;
                else
                    % do nothing because there is no branch between these
                    % two
                end
                
                % and then just index the data from the important slices
                % for each archaeo
                first_outers = outers{i}(top_slice:branch_slice);
                first_inners = inners{i}(top_slice:branch_slice);
                second_outers = outers{j}(top_slice:branch_slice);
                second_inners = inners{j}(top_slice:branch_slice);
                % get those vectors
                [~,first_vec] = centroidline(first_inners,first_outers,0,scale_ratio,0);
                [~,second_vec] = centroidline(second_inners,second_outers,0,scale_ratio,0);
                
                % calculate those angles
                branching_angles(i,j) = branchangle(first_vec,second_vec);
                
                % finally, store those important branching points in 3D for
                % each archaeo (center of archao at branch point
                ellipse = fit_ellipse(first_outers{end}(:,1),first_outers{end}(:,2));
                if isempty(ellipse) || strcmp(ellipse.status, 'Hyperbola found')
                    pgon=polyshape(first_outers{end}(:,1),first_outers{end}(:,2));
                    [x,y]=centroid(pgon);
                else
                    x = ellipse.X0_in;
                    y = ellipse.Y0_in;
                end
                branching_points_3d(i,j,1) = x;
                branching_points_3d(i,j,2) = -y;
                branching_points_3d(i,j,3) = branch_slice * -scale_ratio;
                
            else
                % do nothing
            end
        end
    end
end
