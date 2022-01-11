function [superpix,labeled_pix] = superpixel_trainer(input_ims,n_classes,pix_per_super)
% This function takes a set of input images and does an overpixel
% segmentation on them based upon the approximate number of pixels per
% superpixel input by the user. The function then goes into ginput, and the
% user can click super pixels thatbelong to the number of classes input by
% the user.
%
% IN:
% 
% input_ims: struct where each field contains an image (rgb or grayscale).
% For example input_ims.im1, input_ims.im2, input_ims.im3 could all be
% images. The function will go into each field.
%
% n_classes: The number of classes you would like to gather training data
% for.
%
% pix_per_super: The approximate number of pixels you would like in each
% super pixel to generally define their size.
%
% OUT: 
%
% superpix: struct containing the same fields as input_ims where the value
% of each field is the label matrix given by the superpixels algorithm.
%
% labeled_pix: struct containing the same fields as input_ims where the
% value of each field is a n_class element cell array. Each cell contains a
% column vector containing the indices of the superpixels from superpix
% that correspond to the pixels identified during training.
%
% R. A. Manzuk 11/30/2021
    %% begin the function
    % we know input ims ar coming in as a stuct, good to get the field
    % names
    fields = fieldnames(input_ims);

    % set up an empty struct for the superpixel label mat
    superpix = struct;
    for i = 1:numel(fields)
        % need to know the number of superpixels which is the total number
        % of pixels divided by the pix_per_super
        n_supers = round((size(input_ims.(fields{i}),1) * size(input_ims.(fields{i}),2))/pix_per_super);

        % make the super pixels
        [label_mat,~] = superpixels(input_ims.(fields{i}),n_supers);
        
        % make a boundary mask
        boundaries = boundarymask(label_mat);
        
        % make empty cell array to hold pixel numbers for each class
        class_labels = cell(1,n_classes);
        % now we should ginput the superpixels that represent our classes
        % that we want
        % gotta loop through all classes
        for j = 1:n_classes
            
            % make an empty array to hold the pixel numbers for this class
            this_class = [];

            % show the image 
            imshow(imoverlay(input_ims.(fields{i}),boundaries,'cyan'));
            title("Click superpixels for class " + (j) + ". Press enter to end.");
            hold on
            
            % start a counter
            n=0;
            % make a breakable loop so we can stop clicking
            while true
                % ginput a super pixel
                [col,row] = ginput(1);

                % if we didn't click, break
                if isempty(col) ; break; end

                % move the counter
                n = n+1;

                % figure out which number superpixel was clicked
                clicked_pixel = label_mat(round(row),round(col));

                % update pixel list
                this_class = [this_class;clicked_pixel];

                % we want to indicate this pixel was selected, so need to
                % find its outline and plot a patch over it
                pixel_logical = label_mat == clicked_pixel;
                pixel_boudary = bwboundaries(pixel_logical);
                patch(pixel_boudary{1}(:,2),pixel_boudary{1}(:,1),'white','FaceAlpha',.3);
                drawnow
            end
            %so we can clear the figure for the next class
            hold off

            % put the numbers for this class into the cell array for output
            class_labels{j} = unique(this_class);
        end

        % place the superpixel label mat into the structure for output
        superpix.(fields{i}) = label_mat;
        labeled_pix.(fields{i}) = class_labels;
    end
end