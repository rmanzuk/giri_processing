function [] = plot_class_rectangles(input_im, rectangle_coords, colors)
% this function takes the coordinates of rectangles extracted and plots
% them over the image for visualization.
%
% IN:
% 
% input_im: mxnxn_channel matrix with the image from which the user would
% like the pixel values.
%
% rectangle_coords: matrix that is n_rectangles x 5 where the
% columns are (x_min,x_max,y_min,y_max,class_label). Class label should be
% an integer in the range 1:n_classes. This is the format that the function
% outputs the rectangle coordinates.
%
% colors: (optional) n_classes x 3 matrix containing the RGB color values
% (scaled 0-1) to be used for the boxes for each class plot. If not given,
% the default color order will define the box colors.
%
% OUT: 
% 
% Ryan A. Manzuk 11/16/2021
    %% begin the function   
    % use the final column to tell how many classes there are
    class_list = unique(rectangle_coords(:,end))';

    % did we get colors input? if not, set up colors with the default color
    % order
    if nargin < 3
        color_order = get(0, 'DefaultAxesColorOrder');
        colors = color_order(1:max(class_list),:);
    end

    % show the image
    imshow(input_im)
    hold on
    % iterate through class list and plot the rectangles
    for i = class_list
        % extract just this class' rectangles from the coordinates matrix
        class_inds = rectangle_coords(:,end) == i;
        these_recs = rectangle_coords(class_inds,1:end-1);
        % iterate through and plot each rectangle
        for j = 1:size(these_recs,1)
            % might be easiest to just lay out the points
            x_min = these_recs(j,1);
            x_max = these_recs(j,2);
            y_min = these_recs(j,3);
            y_max = these_recs(j,4);
            plot([x_min, x_min, x_max, x_max, x_min], [y_min, y_max, y_max, y_min, y_min],'Color',colors(i,:))
        end
    end
end