function [pixel_vals,rectangle_coords] = extract_image_rectangles(input_im, n_classes, existing_coords)
% this function goes into an image and extracts the pixel values of
% rectangles for a desired number of classes either through ginput or
% through the input of existing, known rectangle coordinates.
%
% IN:
% 
% input_im: mxnxn_channel matrix with the image from which the user would
% like the pixel values.
%
% n_classes: integer with the number of classes the user would like to
% dessignate.
%
% existing_coords: (optional) matrix that is n_rectangles x 5 where the
% columns are (x_min,x_max,y_min,y_max,class_label). Class label should be
% an integer in the range 1:n_classes. This is the format that the function
% outputs the rectangle coordinates.
%
% OUT: 
%
% pixel_vals: n_pixel x n_channel + 1 matrix with the pixel values for a
% given class (in order of clicking). final column is the class label for
% each pixel
%
% rectangle_coords: matrix that is n_rectangles x 5 where the columns are
% (x_min,x_max,y_min,y_max,class_label). Class label is an integer in the
% range 1:n_classes, given in the order of clicking.
% 
% Ryan A. Manzuk 11/15/2021
    %% begin the function
    % set up emty cell arrays for pixel values and rectangle coordinates
    pixel_vals = [];
    rectangle_coords = [];
    % loop through all classes and get the data via ginput
    for i = 1:n_classes

        % check if the user gave us coordinates or not
        if nargin == 2
            % tell the user what to do
            to_display = 'Please click the rectangle coordinates for class %u\n';
            fprintf(to_display,i);
            
            % if this is the first class, remind user upper left, bottom right
            if i == 1
                fprintf('(click upper left followed by bottom right of each desired rectangle.\n Press enter when finished with class.)\n');
            end

            % show the image and get the user input
            imshow(input_im)
            these_coords = round(ginput);
    
            % rearrange coordinates so the matrix is n_rectangles x 4 and each
            % row contains (x_min,x_max,y_min,y_max)
            rearranged_coords = zeros(size(these_coords,1)/2,4);
            rearranged_coords(:,1) = these_coords(1:2:end-1,1);
            rearranged_coords(:,2) = these_coords(2:2:end,1);
            rearranged_coords(:,3) = these_coords(1:2:end-1,2);
            rearranged_coords(:,4) = these_coords(2:2:end,2);
            % add inclass label column
            rearranged_coords(:,5) = i;
        else % maybe we already got coords and we just have to use them
            class_inds = existing_coords(:,end) == i;
            rearranged_coords = existing_coords(class_inds,:);
        end

        % place in final cell array for rectangle coords
        rectangle_coords = [rectangle_coords;rearranged_coords];

        % and loop through each rectangle and grab the pixels
        pixels = [];
        for j = 1:size(rearranged_coords,1)
            % extract the image rectangle
            image_rec = input_im(rearranged_coords(j,3):rearranged_coords(j,4),rearranged_coords(j,1):rearranged_coords(j,2),:);
            init_pixels = reshape(image_rec,size(image_rec,1)*size(image_rec,2),size(image_rec,3));
            % add the class label column
            class_lables = ones(size(init_pixels,1),1) * i;
            pixels = [pixels;[init_pixels, class_lables]];
        end
        % place in cell array for final export
        pixel_vals = [pixel_vals;pixels];
    end
end