function [inner_3d_rotated, outer_3d_rotated] = rotate_clicked_data3d(inners, outers, block_top_sd, strike_im_heading, bedding_sd,plt)
% This function takes inner and outer clicked data for archaeo branches
% from a GIRI image stack (could be used for other 3D data), along with
% orientation information about the block and field site to rotate the data
% into paleo-gravitational orientation. 
%
% IN
% inners: 1xn_archaeos cell array containing the densified or non-densified 3d outputs for
% inner clicked data from the densify_3d or make_clicking_3d functions.
%
% outers: 1xn_archaeos cell array containing the densified or non-densified 3d outputs for
% outer clicked data from the densify_3d or make_clicking_3d functions.
%
% block_top_sd: 1x2 array containing stike-dip orientation information for
% the top surface of the block or image plane [strike, dip].
%
% strike_im_heading: the north-south-east-west angle heading in degrees of the block
% top arrow of strike in the image plane. Pointing directly to the top of
% the image means an angle of 0 degrees with angles increasing clockwise
% from there.
%
% bedding_sd: 1x2 array containing strike-diop orientation information from the field site of
% the plane of paleo-horizontal (bedding).... [strike, dip].
%
% plt: logical flag if the user would like the 3D rotated data plotted at
% the completion of rotation. 1 for plot, 0 for don't plot
% 
% OUT
% inner_3d_rotated: 1xn_archaeo cell array where each cell contains the 3D
% coordinates for all inner clicked points for a given archaeo. The fourth
% column in the cell is the button data, which may be useful later
%
% outer_3d_rotated: 1xn_archaeo cell array where each cell contains the 3D
% coordinates for all outer clicked points for a given archaeo. The fourth
% column in the cell is the button data, which may be useful later
%
% R. A. Manzuk, 07/24/2020
    %% begin the function
    % we will need several rotation matrices for this....
    % first, rotation matrix for spin around z axis such that strike is facing the top of image
    z_rot_mat = [cosd(strike_im_heading), -sind(strike_im_heading), 0;
        sind(strike_im_heading), cosd(strike_im_heading), 0; 0, 0, 1];

    % then rotation matrix to account for dip (y axis) around
    y_rot_mat = [cosd(-block_top_sd(2)), 0, sind(-block_top_sd(2));
        0, 1, 0; -sind(-block_top_sd(2)), 0, cosd(-block_top_sd(2))];

    % and rotation matrix again around z axis s.t. bedding strike is properly oriented
    net_strike = bedding_sd(1) - block_top_sd(1);
    z_rot_mat2 = [cosd(net_strike), -sind(net_strike), 0;
        sind(net_strike), cosd(net_strike), 0; 0, 0, 1];

    % last rotation matrix....to undo the bedding dip around y axis
    y_rot_mat2 = [cosd(bedding_sd(2)), 0, sind(bedding_sd(2));
        0, 1, 0; -sind(bedding_sd(2)), 0, cosd(bedding_sd(2))];

    % go through, make slice data 3d, and rotate it
    inner_3d_rotated = {};
    for i = 1:numel(inners)
        % extract 3d dense data for this archaeo
        current_points = inners{i};
        % and do all the rotations
        rotated1 = z_rot_mat * [current_points(:,1),current_points(:,2),current_points(:,3)]';
        rotated2 = y_rot_mat * rotated1;
        rotated3 = z_rot_mat2 * rotated2;
        rotated4 = y_rot_mat2 * rotated3;
        final_mat = [rotated4',current_points(:,4)];
        inner_3d_rotated{i} = final_mat;
    end

    % same stuff for the outers data
    outer_3d_rotated = {};
    for i = 1:numel(outers)
        current_points = outers{i};
        % all the rotations
        rotated1 = z_rot_mat * [current_points(:,1),current_points(:,2),current_points(:,3)]';
        rotated2 = y_rot_mat * rotated1;
        rotated3 = z_rot_mat2 * rotated2;
        rotated4 = y_rot_mat2 * rotated3;
        final_mat = [rotated4',current_points(:,4)];
        outer_3d_rotated{i} = final_mat;
    end

if plt
    for i = 1:numel(outer_3d_rotated)
        scatter3(outer_3d_rotated{i}(:,1),outer_3d_rotated{i}(:,2),outer_3d_rotated{i}(:,3),5,'b','filled')
        hold on
        scatter3(inner_3d_rotated{i}(:,1),inner_3d_rotated{i}(:,2),inner_3d_rotated{i}(:,3),5,'r','filled')
    end
else
    % do nothing
end

end




