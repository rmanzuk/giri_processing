function [branch_points_rotated] = rotate_branch_points(branch_points_3d, block_top_sd, strike_im_heading, bedding_sd)
% This function takes inner and outer clicked data for archaeo branches
% from a GIRI image stack (could be used for other 3D data), along with
% orientation information about the block and field site to rotate the data
% into paleo-gravitational orientation. 
%
% IN
% branch_points_3d: n_archaeos x n_archaeos x 3 matrix where the index of
% (i,j,:) gives the 3d location of the center of the ith archaeo at the
% point of branching to the jth archaeo.
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
% branch_points_3d: Just like the input branch points 3d, but rotated to
% account for block orientation and bedding strike/dip.
%
% R. A. Manzuk, 09/18/2020
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
    
    % reshape the branch points matrix to be 3 columns
    x = reshape(branch_points_3d(:,:,1),1,size(branch_points_3d,1)*size(branch_points_3d,2));
    y = reshape(branch_points_3d(:,:,2),1,size(branch_points_3d,1)*size(branch_points_3d,2));
    z = reshape(branch_points_3d(:,:,3),1,size(branch_points_3d,1)*size(branch_points_3d,2));
    % and do all the rotations
    rotated1 = z_rot_mat * [x;y;z];
    rotated2 = y_rot_mat * rotated1;
    rotated3 = z_rot_mat2 * rotated2;
    rotated4 = y_rot_mat2 * rotated3;
    %reshape back to the square matrix style
    x_rotated = reshape(rotated4(1,:),size(branch_points_3d,1),size(branch_points_3d,2));
    y_rotated = reshape(rotated4(2,:),size(branch_points_3d,1),size(branch_points_3d,2));
    z_rotated = reshape(rotated4(3,:),size(branch_points_3d,1),size(branch_points_3d,2));
    % fill in the final mat
    branch_points_rotated = zeros(size(branch_points_3d,1),size(branch_points_3d,2),3);
    branch_points_rotated(:,:,1) = x_rotated;
    branch_points_rotated(:,:,2) = y_rotated;
    branch_points_rotated(:,:,3) = z_rotated;


end




