function [inner_3d_dense,outer_3d_dense] = densify_3d(inners,outers,scale_ratio,densification,plt)
% This function takes inner and outer clicked data for archaeo branches
% from a GIRI image stack (could be used for other 3D data), and
% interpolates between slices to add more points to the point cloud. It
% accomplishes this task by taking the nearest neighbor points from two
% adjacent slices, draws a line between them, and add the desired number of
% points along that line. 

% IN
% inners: 1xn_archaeos cell array containing the cell arrays for inner circles 
% created during data collection. To set this variable up, I just call 
% inners = {inner1, inner2, ..., innern}
%
% outers: 1xn_archaeos cell array containing the cell arrays for outer circles 
% created during data collection. To set this variable up, I just call 
% outers = {outer1, outer2, ..., outern}
%
% scale_ratio: ratio of vertical image separation to pixel-width. For
% example if images are separated by 100 microns and pixels are 20 microns,
% this input is 5.
%
% densification: Factor by which you would like to densify the cloud. For
% example, if you would like 4x as many points, this should be 4, and the
% function will add 3 points for each neighbor pairing between slices.
%
% plt: logical flag if the user would like the 3D rotated data plotted at
% the completion of rotation. 1 for plot, 0 for don't plot
% 
% OUT
% inner_3d_dense: 1xn_archaeo cell array where each cell contains the 3D
% coordinates for all inner clicked points for a given archaeo plus the added points for density. The fourth
% column in the cell is the button data, which may be useful later
%
% outer_3d_rotated: 1xn_archaeo cell array where each cell contains the 3D
% coordinates for all outer clicked points for a given archaeo plus the added points for density. The fourth
% column in the cell is the button data, which may be useful later
%
% R. A. Manzuk, 08/11/2020
    %% begin the function
    % turn the densification factor into number of points that need to be
    % added between slices
    num_points = densification -1;
    
    % go through each archaeo's inners cell array and densify
    inner_3d_dense = {};
    for i = 1:numel(inners)
        current_points = [];
        % need to go through slice by slice
        for j = 1:numel(inners{i})-1
            % define this slice and the next in 3d space given the spacing
            % between slices
            this_slice = [inners{i}{j},ones(size(inners{i}{j},1),1).*(j * -scale_ratio)];
            next_slice = [inners{i}{j+1},ones(size(inners{i}{j+1},1),1).*((j+1) * -scale_ratio)];
            % if the current slice isn't empty, might as well extract its
            % points
            if ~isempty(this_slice)
                current_points = [current_points; this_slice];
            else
                %do nothing
            end
            % then if both slices contain points, look to interpolate
            % between them
            if ~isempty(this_slice) && ~isempty(next_slice)
                % find the point in the next slice that is closest to each
                % point in the current slice
                nearest = dsearchn([next_slice(:,1),next_slice(:,2),next_slice(:,4)], [this_slice(:,1),this_slice(:,2),this_slice(:,4)]);
                % go through each point match, and add the desired number
                % of points along a line connecting the two. this is our
                % interpolation.
                for k = 1:length(nearest)
                    point1 = this_slice(k,[1,2,4]);
                    point2 = next_slice(nearest(k),[1,2,4]);
                    x_step = (point2(1) - point1(1))/(densification+1); 
                    y_step = (point2(2) - point1(2))/(densification+1); 
                    z_step = (point2(3) - point1(3))/(densification+1); 
                    new_points = [];
                    for l = 1:num_points
                        new_points(l,:) = [point1(1) + (l*x_step),point1(2) + (l*y_step),point1(3) + (l*z_step)];
                    end
                    % need to account for button data
                    fake_button = zeros(size(new_points,1),1);
                    final_points = [new_points(:,1),new_points(:,2),fake_button,new_points(:,3)];
                    current_points = [current_points; final_points];
                end
            else
                %do nothing
            end
        end
        % the last slice just needs to be added on
        last_slice = [inners{i}{end},ones(size(inners{i}{end},1),1).*(numel(inners{i}) * -scale_ratio)];
        current_points = [current_points; this_slice];
        % then just reorder so button data is last column, and account for
        % ginput y -1
        inner_3d_dense{i} = [current_points(:,1),current_points(:,2).*-1,current_points(:,4),current_points(:,3)];
    end

    % and do the same thing for outer data
    outer_3d_dense = {};
    for i = 1:numel(outers)
        current_points = [];
        for j = 1:numel(outers{i})-1
            this_slice = [outers{i}{j},ones(size(outers{i}{j},1),1).*(j * -scale_ratio)];
            next_slice = [outers{i}{j+1},ones(size(outers{i}{j+1},1),1).*((j+1) * -scale_ratio)];
            if ~isempty(this_slice)
                current_points = [current_points; this_slice];
            else
                %do nothing
            end
            if ~isempty(this_slice) && ~isempty(next_slice)
                nearest = dsearchn([next_slice(:,1),next_slice(:,2),next_slice(:,4)], [this_slice(:,1),this_slice(:,2),this_slice(:,4)]); 
                for k = 1:length(nearest)
                    point1 = this_slice(k,[1,2,4]);
                    point2 = next_slice(nearest(k),[1,2,4]);
                    x_step = (point2(1) - point1(1))/(densification+1); 
                    y_step = (point2(2) - point1(2))/(densification+1); 
                    z_step = (point2(3) - point1(3))/(densification+1); 
                    new_points = [];
                    for l = 1:num_points
                        new_points(l,:) = [point1(1) + (l*x_step),point1(2) + (l*y_step),point1(3) + (l*z_step)];
                    end
                    fake_button = zeros(size(new_points,1),1);
                    final_points = [new_points(:,1),new_points(:,2),fake_button,new_points(:,3)];
                    current_points = [current_points; final_points];
                end
            else
                %do nothing
            end
        end
        last_slice = [outers{i}{end},ones(size(outers{i}{end},1),1).*(numel(outers{i}) * -scale_ratio)];
        current_points = [current_points; this_slice];
        outer_3d_dense{i} = [current_points(:,1),current_points(:,2).*-1,current_points(:,4),current_points(:,3)];
    end

    % plot the whole thing if we want
    if plt
        for i = 1:numel(inner_3d_dense)
            scatter3(outer_3d_dense{i}(:,1),outer_3d_dense{i}(:,2),outer_3d_dense{i}(:,3),5,'b','filled')
            hold on
            scatter3(inner_3d_dense{i}(:,1),inner_3d_dense{i}(:,2),inner_3d_dense{i}(:,3),5,'r','filled')
            hold on
        end
    end

end