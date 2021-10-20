function [sharpness_vector, upper_lefts, window_row_top] = sharpness_dissipation_curve(img, window_size, window_stride)
% This function takes an image and calculates the sharpness moving from the
% center to the right edge for windows of a defined stride and size. The
% goal is to quantify the sharpness dissipation curve of the lens.
%
% IN
%
% img: 2D matrix containing the image in any format.
%
% window_size: single value representing size of the windows in pixels. 
% Windows will be square.
%
% window stride: single value representing number of pixels the window is
% moved between each sampling
%
% OUT
%
% sharpness_vector: 1xn_windows vector containing the sharpness values of
% the windows moving from center to the edge.
%
% upper_lefts: 1xn_windows vector containing column coordinates for the 
% left edge of each window moving from center to edge.
%
% window_row_top: row index of the top of all windows
%
% 06/10/2021 Ryan A. Manzuk
    %% begin the function
    % figure out the column index where the first window should start
    window_col_start = (size(img,2)/2) - window_size/2;

    % and the row index that is the top of all boxes
    window_row_top = (size(img,1)/2) - window_size/2;

    % where is the last box we can make without going off the image
    window_max = size(img,2) - window_size;

    % use that information to make a vector of column indices for the starts of
    % each box based upon the stride
    upper_lefts = [window_col_start:window_stride:window_max];

    % calculate the image gradient
    [grad_x,grad_y] = gradient(img);
    grad_norms = sqrt(grad_x.*grad_x+grad_y.*grad_y);

    % set up the sharpness vector
    sharpness_vector = zeros(1,numel(upper_lefts));

    % iterate through the windows and get the sharpnesses
    for i = 1:numel(upper_lefts)
        this_grad = grad_norms(window_row_top:(window_row_top+window_size), upper_lefts(i):(upper_lefts(i)+window_size));
        sharpness_vector(i) = sum(sum(this_grad))./numel(this_grad);
    end

end