function [opened_classified] = multiclass_open(classified_im,rad,n_pix)
% This function takes a fully classified, indexed image and performs a
% morphological opening on all classes given a certain size of island (in
% number of pixels) to be removed. The islands are filled by the
% surrounding class.
%
% IN:
% 
% classified_im: indexed classified image (mxn integer matrix).
%
% rad: radius of the opening and closing disc
%
% OUT:
%
% opened_classified: indexed classified image with islands removed, filled
% by the class that surrounds them.
%
% R. A. Manzuk 11/30/2021
%
    %% begin the function
    % make an empty image to hold the final classified, opened image
    opened_classified = zeros(size(classified_im));

    % define the opening/closing disc
    SE = strel('disk',rad);

    % loop through all of the classes
    for i = 1:max(unique(classified_im))
        % grab just a logical of this class
        this_logical = classified_im == i;
        
        % open it to get rid of smaller elements of specified size
        opened_logical = imopen(this_logical,SE);

        % and close it so it can swallow the small elements that are
        % removed from other classes
        closed_logical = imclose(opened_logical,SE);

        % end by placing it in the final classification to be output.
        opened_classified(closed_logical) = i;
    end

    % some parts are left as 0, so need to go in and fill them with the
    % pixel value to the right.

    zero_logical = opened_classified == 0;
    % get the connected components of the 0 spots
    concomps = bwconncomp(zero_logical);

    % loop through each and fill them in
    for i = 1:numel(concomps.PixelIdxList)
        % get the indeces of the pixels concerned
        [row,col] = ind2sub(size(zero_logical),concomps.PixelIdxList{i});
        
        % let's just use the median row
        med_row = round(median(row));

        % and figure out which column is immediately to the right of this
        % patch in that row
        next_col = max(col(row==med_row)) + 1;

        % if this patch is all the way to the right, need to use the left
        if max(col(row==med_row)) == size(zero_logical,2)
            next_col = min(col(row==med_row)) - 1;
        end

        % finally fill in 
        opened_classified(row,col) = opened_classified(med_row,next_col);
    end

    classified_im = opened_classified;
    opened_classified = zeros(size(classified_im));

    % loop through all of the classes
    for i = 1:max(unique(classified_im))
        % grab just a logical of this class
        this_logical = classified_im == i;
        
        % open it to get rid of smaller elements of specified size
        opened_logical = bwareaopen(this_logical,n_pix);

        % and close it so it can swallow the small elements that are
        % removed from other classes
        closed_logical = ~bwareaopen(~opened_logical,n_pix);

        % end by placing it in the final classification to be output.
        opened_classified(closed_logical) = i;
    end

    % some parts are left as 0, so need to go in and fill them with the
    % pixel value to the right.

    zero_logical = opened_classified == 0;
    % get the connected components of the 0 spots
    concomps = bwconncomp(zero_logical);

    % loop through each and fill them in
    for i = 1:numel(concomps.PixelIdxList)
        % get the indeces of the pixels concerned
        [row,col] = ind2sub(size(zero_logical),concomps.PixelIdxList{i});
        
        % let's just use the median row
        med_row = round(median(row));

        % and figure out which column is immediately to the right of this
        % patch in that row
        next_col = max(col(row==med_row)) + 1;

        % if this patch is all the way to the right, need to use the left
        if max(col(row==med_row)) == size(zero_logical,2)
            next_col = min(col(row==med_row)) - 1;
        end

        % finally fill in 
        opened_classified(row,col) = opened_classified(med_row,next_col);
    end

end