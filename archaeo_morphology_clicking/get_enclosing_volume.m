function [convhull_points, enclosing_volume] = get_enclosing_volume(outers_3d, unit_conversion)
% This function takes 3d data of objects arranged in cell arrays and
% calculates the 3d volume of their enclosing convex hull
%
% IN
% outers_3d: 1xn_branches cell array containing the densified or non-densified 3d outputs for
% outer clicked data from the densify_3d or make_clicking_3d functions.
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

% R. A. Manzuk, 03/23/2021
    %% begin the function
    
    
    % set up empty array, and just go through the cell array and stack all 
    %points into one 3-column matrix
    all_3d = [];
    for i = 1:numel(outers_3d)
        all_3d = [all_3d; outers_3d{i}(:,1:3).*unit_conversion];
    end
    
    % then the enclosing volume is the convex hull of those 3d points.
    [outer_points,enclosing_volume] = convhull(all_3d);
    convhull_points = [all_3d(outer_points,1),all_3d(outer_points,2)];
        
end