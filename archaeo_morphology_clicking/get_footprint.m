function [footprint_points, footprint_area] = get_footprint(outers_3d, dim, unit_conversion)
% This function takes 3d data of objects arranged in cell arrays and
% calculates the 2d footprint of all combined data on a desired plane
%
% IN
% outers_3d: 1xn_branches cell array containing the densified or non-densified 3d outputs for
% outer clicked data from the densify_3d or make_clicking_3d functions.
%
% dim: dimension that should be looked along to get the footprint. For
% example, if the user wants to look down the z-axis to get the xy plane
% footprint, this input should be 3.
%
% unit_conversion: multiplication factor to take the units from pixels to
% whatever unit you want to measure in. For example, if each pixel
% represents 500 um and you want to measure in cm, this number would be
% 0.05.
% 
% OUT
% footprint_points: n_points x 2 matrix with the 2d points that comprise
% the outline of the footrpint.
%
% footprint_area: area of the footprint in square of the input units

% R. A. Manzuk, 02/25/2021
    %% begin the function
    
    % the dimensions we want in the footprint are the set difference
    % between 1,2,3 and the input dimension.
    dims_take = setdiff([1,2,3],dim);
    
    % set up empty array, and just go through the cell array and stack all 
    %points into one 2d matrix
    all_2d = [];
    for i = 1:numel(outers_3d)
        all_2d = [all_2d; outers_3d{i}(:,dims_take).*unit_conversion];
    end
    
    % then the footprint is the convex hull of those 2d points.
    [outer_points,footprint_area] = convhull(all_2d);
    footprint_points = [all_2d(outer_points,1),all_2d(outer_points,2)];
        
end