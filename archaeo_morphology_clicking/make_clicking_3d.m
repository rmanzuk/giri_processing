function [outlines_3d] = make_clicking_3d(outlines,scale_ratio)
% This function just takes the data from clicking (dense slices or not) and puts it into 3d space
% for later processing. The point cloud is not further densified or altered in any
% way
%
% IN
% outlines: 1xn_branches cell array containing the cell arrays for outer circles 
% created during data collection.
%
% scale_ratio: ratio of vertical image separation to pixel-width. For
% example if images are separated by 100 microns and pixels are 20 microns,
% this input is 5.
%
% 
% OUT
%
% outlines_3d: 1xn_branch cell array where each cell contains the 3D
% coordinates for all outer clicked points for a given archaeo. The fourth
% column in the cell is the button data, which may be useful later
%
% Ryan A. Manzuk 08/12/2020
    %% begin the function
    % go through each archaeo individually
    outlines_3d = {};
    for i = 1:numel(outlines)
        current_points = [];
        % need to go through slice by slice
        for j = 1:numel(outlines{i})
            % define this slice in 3d space given the spacing between slices
            this_slice = [outlines{i}{j},ones(size(outlines{i}{j},1),1).*(j * -scale_ratio)];
            % if the current slice isn't empty, might as well extract its
            % points
            if ~isempty(this_slice)
                current_points = [current_points; this_slice];
            else
                %do nothing
            end
        end
        % then just reorder so button data is last column, and account for
        % ginput y -1
        outlines_3d{i} = [current_points(:,1),current_points(:,2).*-1,current_points(:,4),current_points(:,3)];
    end
    
end