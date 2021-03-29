function [branched_flags,branching_points_3d] = process_branched_network(outlines,scale_ratio)
% This function takes a whole slew of archaeo data and figures out which
% archaeos branch to each other and what the angle of branching is. Right
% now, the the indicator of branching is having overlapping outer circle
% datapoints. That method could easily be changed if another were thought
% up. 
%
% IN 
% outlines: 1xn_branches cell array containing the outlines,
% separated slice-wise, created though automated or clicking data
% collection. Could be densified.
%
% scale_ratio: ratio of vertical image separation to pixel-width. For
% example if images are separated by 100 microns and pixels are 20 microns,
% this input is 5.
%
% OUT
% branched_flags: n_archaeos x n_archaeos comparison matrix where branches
% are kept track of. For any non-zero value in this matrix, the indices
% indicate which archaeos branch, and the values at those indices indicate
% which slice contains the branching point. Because each archaeo is
% considdered in both the rows and columns, this matrix is redundant you
% really only need the half of the data above (or below) the diagonal.
% Think of it like a covariance matrix. 
%
% branch_points_3d: n_archaeos x n_archaeos x 3 matrix where the index of
% (i,j,:) gives the 3d location of the center of the ith archaeo at the
% point of branching to the jth archaeo.

% R. A. Manzuk 07/18/2020

    %% begin dat sweet function
    % how many slices are we working with?
    n_slices = max([cellfun(@numel,outlines)],[],'all');

    % make sure all individual archaeo cell arrays are same size...just in
    % case

    for i = 1:numel(outlines)
        if numel(outlines{i}) < n_slices
            add_empty = cell(1,n_slices - numel(outlines{i}));
            outlines{i}(numel(outlines{i})+1:n_slices) = add_empty;
        else
            % do nothing
        end
    end

    %how many archaeos?
    n_archaeos = numel(outlines);

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
        slice_outers = {};
        for j = 1:n_archaeos
            slice_outers(j) = outlines{j}(i);
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
    % number this happens. Now, we just where these branching points are in 3d.
    
    % Keep track of branch points in a similar comparison matrix where in stead of
    % noting slice number of the branching point, we'll note a 3d position.
    branching_points_3d = zeros(n_archaeos,n_archaeos,3);
    
    % Go through each entry of the flagging matrix
    for i = 1:n_archaeos
        for j = 1:n_archaeos
            % if it's not zero, we know we've got a branch
            if branched_flags(i,j) ~= 0
                % finally, store those important branching points in 3D for
                % each archaeo (center of archao at branch point)
                
                % index the first branch's outline at the branch point
                branch_slice = branched_flags(i,j);
                first_outers = outlines{i}{branch_slice};
                % fit ellipse (or polygon)
                ellipse = fit_ellipse(first_outers(:,1),first_outers(:,2));
                if isempty(ellipse) || strcmp(ellipse.status, 'Hyperbola found')
                    pgon=polyshape(first_outers(:,1),first_outers(:,2));
                    [x,y]=centroid(pgon);
                else
                    x = ellipse.X0_in;
                    y = ellipse.Y0_in;
                end
                % take the center
                branching_points_3d(i,j,1) = x;
                branching_points_3d(i,j,2) = -y;
                branching_points_3d(i,j,3) = branch_slice * -scale_ratio;
                
            else
                % do nothing
            end
        end
    end
end
