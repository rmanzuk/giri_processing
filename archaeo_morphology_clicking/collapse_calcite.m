function [collapsed_outers, collapsed_inners, collapsed_branch_points] = collapse_calcite(crack_left_3d, crack_right_3d, outers_3d, inners_3d, branch_points_3d)
% this function takes in 3d data for branching archaeos and 3d data for any
% crack or calcite that might need to be removed and adjusted for. It
% returns new 3d data of for the 3d locations of teh branched archaeo point
% clouds based upon the removal of the crack. 
%
% IN
% crack_left_3d: 1x1 cell array containing densified or non-densifiec 3d
% outputs the clicked points of the left edge of the crack/calcite.
%
% crack_right_3d: 1x1 cell array containing densified or non-densifiec 3d
% outputs the clicked points of the right edge of the crack/calcite.
%
% outers_3d: 1xn_branches cell array containing the densified or
% non-densified 3d outputs for outer clicked data from the densify_3d or
% make_clicking_3d functions.
%
% inners_3d: 1xn_branches cell array containing the densified or
% non-densified 3d outputs for inner clicked data from the densify_3d or
% make_clicking_3d functions.
%
% branch_points_3d: n_archaeos x n_archaeos x 3 matrix where the index of
% (i,j,:) gives the 3d location of the center of the ith archaeo at the
% point of branching to the jth archaeo. Output from the process branched
% network function.
%
% OUT
% 
% collapsed_outers: 1xn_branches cell array containint the 3d data for each
% branch's outers with the crack/calcite collapsed. Still in the original
% coordinte system of the input data.
%
% collapsed_inners: 1xn_branches cell array containint the 3d data for each
% branch's inners with the crack/calcite collapsed. Still in the original
% coordinte system of the input data.
%
% collapsed_branch_points: same as the input branch_points_3d, but with any
% branch points moved to account for the slight movement of branches when
% collapsing the calcite/crack.
%
% Ryan Manzuk 03/01/2021
%% Begin the function
    % assemble all the crack points into one matrix
    all_crack_points = crack_left_3d{1}(:,1:3);
    all_crack_points = [all_crack_points; crack_right_3d{1}(:,1:3)];

    % regress a single plane over the crack
    B = [all_crack_points(:,1),all_crack_points(:,2),ones(size(all_crack_points,1),1)] \ all_crack_points(:,3);
    
    % from the equation of the plane, figure out the inclination and
    % declination of the normal vector
    [declination,slope_run] = cart2pol(B(1),B(2));
    inclination = atand(1/slope_run);
    declination = rad2deg(declination);
    declination = -1.*(declination - 90);

    if declination < 0
        declination = 360 + declination;
    end

    % based upon the normal inclination and declination, set up rotation
    % matrices that can be used to rotate the crack and branches into a
    % space where the calcite plane is horizontal
    z_rot_mat = [cosd(declination), -sind(declination), 0;
                sind(declination), cosd(declination), 0; 0, 0, 1];
    x_rot_mat = [1,0,0;0,cosd(90-inclination), sind(90-inclination);
                0, -sind(90-inclination), cosd(90-inclination)];
   
    % rotate the left crack
    left_crack_rot1 = z_rot_mat * crack_left_3d{1}(:,1:3)';
    left_crack_rotated = x_rot_mat * left_crack_rot1;
    left_crack_rotated = left_crack_rotated';

    % rotate the right crack
    right_crack_rot1 = z_rot_mat * crack_right_3d{1}(:,1:3)';
    right_crack_rotated = x_rot_mat * right_crack_rot1;
    right_crack_rotated = right_crack_rotated';
    
    % now, in the rotated space of the crack we'll want the top and bottom
    % represented by gridded surfaces. Start by setting up a meshgrid over
    % the x,y coords occupied by the crack.
    rotated_crack_all = [left_crack_rotated; right_crack_rotated];
    minima = min(rotated_crack_all);
    maxima = max(rotated_crack_all);
    warning('off','all')
    [X,Y] = meshgrid([minima(1):10:maxima(1)],[minima(2):10:maxima(2)]);
       
    % then just grid the crack surface data in the X Y coordinate space
    gridded_left = griddata(left_crack_rotated(:,1),left_crack_rotated(:,2),left_crack_rotated(:,3),X,Y);
    gridded_right = griddata(right_crack_rotated(:,1),right_crack_rotated(:,2),right_crack_rotated(:,3),X,Y);
    warning('on','all')
    
    % now go through the branches and rotate them into the space with the
    % horizontal crack
    collapsed_outers = {};
    collapsed_inners = {};
    collapsed_branch_points = zeros(size(branch_points_3d));
    for i = 1:numel(outers_3d)
        
        % grab all of the necessary components
        these_outers = outers_3d{i}(:,1:3);
        these_inners = inners_3d{i}(:,1:3); 
        these_branch_points = squeeze(branch_points_3d(i,:,:));
        these_branch_points = these_branch_points(any(these_branch_points,2),:);

        % rotate everything s.t. they are in the same space as the
        % horizontal crack.
        outers_rot1 = z_rot_mat * these_outers';
        outers_rotated = x_rot_mat * outers_rot1;
        outers_rotated = outers_rotated';

        inners_rot1 = z_rot_mat * these_inners';
        inners_rotated = x_rot_mat * inners_rot1;
        inners_rotated = inners_rotated';
        
        % only rotate branch points if they exist.
        if ~isempty(these_branch_points)
            branch_points_rot1 = z_rot_mat * these_branch_points';
            branch_points_rotated = x_rot_mat * branch_points_rot1;
            branch_points_rotated = branch_points_rotated';
        else
            branch_points_rotated = [];
        end

        % now figure out which part of the calcite we need to consider for this
        % branch i.e. what part is the branch closest to the calcite

        [minima,~] = min(pdist2(right_crack_rotated,outers_rotated));
        [~,closest] = min(minima);

        % based upon that closest point's z position, and a range of 20 above
        % it, figure out the range of xy values to consider (this will be the
        % local area of the crack we check.
        z_range = [outers_rotated(closest,3),outers_rotated(closest,3)+20];
        points_in_zrange = outers_rotated(:,3) > z_range(1) & outers_rotated(:,3) < z_range(2);

        % given the branch points in the z range, what is the range of x,y
        % values we should consider
        min_x = min(outers_rotated(points_in_zrange,1));
        max_x = max(outers_rotated(points_in_zrange,1));
        min_y = min(outers_rotated(points_in_zrange,2));
        max_y = max(outers_rotated(points_in_zrange,2));

        % and use those ranges to make a logical of the crack's xy grid
        if isempty(max_x)
            mean_calcite_upper = mean(gridded_right,'all','omitnan');
            mean_calcite_lower = mean(gridded_left,'all','omitnan');
        else
            crack_y_inds = find(X(1,:) < max_x & X(1,:) > min_x); 
            crack_x_inds = find(Y(:,1) < max_y & Y(:,1) > min_y);

            % find the mean upper and lower height of the calcite at the positions
            % indicated by this branch
            mean_calcite_upper = mean(gridded_right(crack_x_inds,crack_y_inds),'all','omitnan');
            mean_calcite_lower = mean(gridded_left(crack_x_inds,crack_y_inds),'all','omitnan');
            
            if isnan(mean_calcite_upper) || isnan(mean_calcite_upper)
                mean_calcite_upper = mean(gridded_right,'all','omitnan');
                mean_calcite_lower = mean(gridded_left,'all','omitnan');
            end
  
        end

        % then define points above, within, and below the calcite
        points_above = outers_rotated(:,3) >= mean_calcite_upper;
        points_below = outers_rotated(:,3) <= mean_calcite_lower;

        % assess any offset created by the crack, just through the position of
        % the points immediately above and below the crack.
        just_above = points_above & outers_rotated(:,3) < mean_calcite_upper + 50;
        just_below = points_below & outers_rotated(:,3) > mean_calcite_lower - 100;

        mean_xy_above = mean(outers_rotated(just_above,1:2));
        mean_xy_below = mean(outers_rotated(just_below,1:2));
        
        if isnan(sum(mean_xy_above)) || isnan(sum(mean_xy_below))
            mean_xy_above = [0,0];
            mean_xy_below = [0,0];
        end

        below_translation = [2*(mean_xy_above - mean_xy_below)/3,mean_calcite_upper-mean_calcite_lower];

        new_outers = [outers_rotated(points_above,:);(outers_rotated(points_below,:) + below_translation)];

        % and apply translation to inners and branch points
        inners_above = inners_rotated(:,3) >= mean_calcite_upper;
        inners_below = inners_rotated(:,3) <= mean_calcite_lower;
        new_inners = [inners_rotated(inners_above,:);(inners_rotated(inners_below,:) + below_translation)];
        
        % only worry about branch points if they're there.
        if ~isempty(branch_points_rotated)
            bps_above = branch_points_rotated(:,3) >= mean_calcite_upper;
            bps_below = branch_points_rotated(:,3) < mean_calcite_upper;

            new_bps = [branch_points_rotated(bps_above,:);(branch_points_rotated(bps_below,:) + below_translation)];
        end

        % rotate everybody back to the original orientation
        z_rot_mat2 = [cosd(declination), sind(declination), 0;
                -sind(declination), cosd(declination), 0; 0, 0, 1];
        x_rot_mat2 = [1,0,0;0,cosd(90-inclination), -sind(90-inclination);
                0, sind(90-inclination), cosd(90-inclination)];

        new_out_rot1 = x_rot_mat2 * new_outers';
        final_outers = z_rot_mat2 * new_out_rot1;
        collapsed_outers{i} = final_outers';

        new_inn_rot1 = x_rot_mat2 * new_inners';
        final_inners = z_rot_mat2 * new_inn_rot1;
        collapsed_inners{i} = final_inners';
        
        % handle branch points and their messy indexing
        if ~isempty(branch_points_rotated)
            new_bps_rot1 = x_rot_mat2 * new_bps';
            final_bps = z_rot_mat2 * new_bps_rot1;
            final_bps = final_bps';

            indices = find(branch_points_3d(i,:,1) ~= 0);
            for j = 1:length(indices)
                collapsed_branch_points(i,indices(j),:) = reshape(final_bps(j,:),1,1,3);
            end
        end
    end
end

