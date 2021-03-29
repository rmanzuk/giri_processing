function [data_rotated] = rotate_clicked_data3d(data_3d, block_top_sd, strike_im_heading, bedding_sd)
% This function takes inner and outer clicked data for archaeo branches
% from a GIRI image stack (could be used for other 3D data), along with
% orientation information about the block and field site to rotate the data
% into paleo-gravitational orientation. 
%
% IN
% data_3d: 1xn_branch cell array 3d point data of each branch. Could be
% center lines or outlines.
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
% OUT
% data_rotated: 1xn_branch cell array where each cell contains the 3D
% coordinates for all input data rotated based upon input orientations. The fourth
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
    data_rotated = {};
    for i = 1:numel(data_3d)
        if ~isempty(data_3d{i})
            % extract 3d dense data for this archaeo
            current_points = data_3d{i};
            % and do all the rotations
            rotated1 = z_rot_mat * [current_points(:,1),current_points(:,2),current_points(:,3)]';
            rotated2 = y_rot_mat * rotated1;
            rotated3 = z_rot_mat2 * rotated2;
            rotated4 = y_rot_mat2 * rotated3;
            if size(current_points,2) == 4
                final_mat = [rotated4',current_points(:,4)];
            else
                final_mat = rotated4';
            end
            data_rotated{i} = final_mat;
        else
            data_rotated{i} = [];
        end
    end
end




