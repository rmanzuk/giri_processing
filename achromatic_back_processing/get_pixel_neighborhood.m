function [neighborhoods] = get_pixel_neighborhoods(image,grid_coords,n)
% function to extract a neighboorhood of 2n+1 pixels centered around the
% coornates specified for one or multiple points

% Note;

% IN:
% image: the image from which you want to extract neighborhoods
% grid_coords: 2-column matrix containing the row and column pixel
% positions where neighborhoods should be extracted from the image
% n: size of the neighborhood, where height and width will be 2n + 1,
% centered around the pixel specified in grid_coords

% OUT;
% neighborhoods: 3-channel image neighborhoods of specified size and
% position

% R. Manzuk
% 23 January, 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % make sure that image is a double
    image = im2double(image);
    
    % define point to upper left and lower right of pixel coordinate to create
    % our definition of the neighborhood
    lower_rights = [grid_coords(:,1) - n,grid_coords(:,2) + n];
    upper_lefts = [grid_coords(:,1) + n,grid_coords(:,2) - n];

    %then extract
    neighborhoods = [];
    for i = 1:size(upper_lefts,1)
        neighborhoods(:,:,:,i) = image((round(lower_rights(i,1)):round(upper_lefts(i,1))),(round(upper_lefts(i,2)):round(lower_rights(i,2))),:); 
    end
end