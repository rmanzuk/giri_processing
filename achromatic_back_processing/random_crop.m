function[crops, crop_coords] = random_crop(input_im,n_crops,crop_size)
% This function goes through and makes a desired number of random crops of
% an image of a certain size and puts them into a struct.
%
% IN:
%
% input_im: mxnxn_channel matrix with the image from which the user would
% like the crops.
%
% n_crops: how many crops the user would like.
%
% 1x2 vector containing the row and column dimensions of the
% desired crops to be sampled.
%
% OUT:
%
% crops: struct containing the crops in the fields 'crop1',
% 'crop2',...,'cropn'.
%
% crop_coords: n_class x 4 matrix with the coordiates for the crops on the
% image where each row is in the format: x_min, x_max, y_min, y_max.
%
% 12/01/2021 R.A. Manzuk
    %% begin the function
    % get the randomly placed upper corners
    crop_corners = round(rand(n_crops,2) .* (size(input_im,[1,2]) - crop_size-1));

    % make an empty struct for the crops
    crops = struct;
    
    % empty array to grab the crop coordinates
    crop_coords = zeros(n_crops,4);
    % we can fill that array right away
    crop_coords(:,1) = crop_corners(:,2);
    crop_coords(:,2) = crop_corners(:,2) + crop_size(2);
    crop_coords(:,3) = crop_corners(:,1);
    crop_coords(:,4) = crop_corners(:,1) + crop_size(1); 
    % make the crops
    for i = 1:n_crops
        % crop it
        this_bit = input_im(crop_corners(i,1):crop_corners(i,1)+crop_size(1)-1,...
        crop_corners(i,2):crop_corners(i,2)+crop_size(2)-1,:);
        
        % make the field for the struct
        field = strcat('crop',num2str(i));

        % put it in the struct
        crops.(field) = this_bit;
    end
end