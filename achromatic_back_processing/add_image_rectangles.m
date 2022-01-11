function [pixel_vals,rectangle_coords] = add_image_rectangles(input_im,init_pixel_vals,init_rectangle_coords,existing_coords)
% this function adds rectangles to previously existing ones when using them
% to gather training data
%
% IN:
% 
% input_im: mxnxn_channel matrix with the image from which the user would
% like the pixel values.
%
% init_pixel_vals: output pixel values from the original selecting of
% rectangles.
%
% init_rectangle_coords: output rectangle coordinates from the original
% selecting of rectangles.
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
% Ryan A. Manzuk 11/29/2021
    %% begin the function
    % set up emty cell arrays for pixel values and rectangle coordinates
    pixel_vals = init_pixel_vals;
    rectangle_coords = init_rectangle_coords;
    n_classes = max(unique(init_rectangle_coords(:,end)));
    % loop through all classes and get the data via ginput
    for i = 1:n_classes
        % grab the indeces of this class in the existing pixel vals and
        % rectangle coords
        class_rec_inds = rectangle_coords(:,end) == i;

        % check if the user gave us coordinates or not
        if nargin == 3
            % tell the user what to do
            to_display = 'Please click the rectangle coordinates for class %u\n';
            fprintf(to_display,i);
            
            % if this is the first class, remind user upper left, bottom right
            if i == 1
                fprintf('(click upper left followed by bottom right of each desired rectangle.\n Press enter when finished with class.)\n');
            end

            % show the image with existing rectangles plotted on top
            plot_class_rectangles(input_im(:,:,1:3),rectangle_coords(class_rec_inds,:));
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

    % just sort at the end to make nice
    rectangle_coords = sortrows(rectangle_coords,size(rectangle_coords,2));
    pixel_vals = sortrows(pixel_vals,size(pixel_vals,2));
end