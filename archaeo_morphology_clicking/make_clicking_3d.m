function [inner_3d,outer_3d] = make_clicking_3d(inners,outers,scale_ratio,plt)
% This function just takes the data from clicking (dense slices or not) and puts it into 3d space
% for later processing. The point cloud is not further densified or altered in any
% way

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
%
% plt: logical flag if the user would like the 3D rotated data plotted at
% the completion of rotation. 1 for plot, 0 for don't plot
% 
% OUT
% inner_3d: 1xn_archaeo cell array where each cell contains the 3D
% coordinates for all inner clicked points for a given archaeo . The fourth
% column in the cell is the button data, which may be useful later
%
% outer_3d_rotated: 1xn_archaeo cell array where each cell contains the 3D
% coordinates for all outer clicked points for a given archaeo. The fourth
% column in the cell is the button data, which may be useful later
%
% Ryan A. Manzuk 08/12/2020
    %% begin the function
    % go through each archaeo individually
    inner_3d = {};
    for i = 1:numel(inners)
        current_points = [];
        % need to go through slice by slice
        for j = 1:numel(inners{i})
            % define this slice in 3d space given the spacing between slices
            this_slice = [inners{i}{j},ones(size(inners{i}{j},1),1).*(j * -scale_ratio)];
            % if the current slice isn't empty, might as well extract its
            % points
            if ~isempty(this_slice)
                current_points = [current_points; this_slice];
            else
                %do nothing
            end
        end
        % then just reorder so button data is last column
        inner_3d{i} = [current_points(:,1),current_points(:,2),current_points(:,4),current_points(:,3)];
    end
    
    % and do the same thing for outers
    outer_3d = {};
    for i = 1:numel(outers)
        current_points = [];
        % need to go through slice by slice
        for j = 1:numel(outers{i})
            % define this slice in 3d space given the spacing between slices
            this_slice = [outers{i}{j},ones(size(outers{i}{j},1),1).*(j * -scale_ratio)];
            % if the current slice isn't empty, might as well extract its
            % points
            if ~isempty(this_slice)
                current_points = [current_points; this_slice];
            else
                %do nothing
            end
        end
        % then just reorder so button data is last column
        outer_3d{i} = [current_points(:,1),current_points(:,2),current_points(:,4),current_points(:,3)];
    end
    
    
    % plot the whole thing if we want
    if plt
        for i = 1:numel(inner_3d)
            scatter3(outer_3d{i}(:,1),outer_3d{i}(:,2),outer_3d{i}(:,3),5,'b','filled')
            hold on
            scatter3(inner_3d{i}(:,1),inner_3d{i}(:,2),inner_3d{i}(:,3),5,'r','filled')
            hold on
        end
    end
end