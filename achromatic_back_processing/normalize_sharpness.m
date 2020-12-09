function [global_normalized,window_normalized] = normalize_sharpness(global_sharpness,window_sharpness)
% This function takes the output sharpness metrics for images and windows
% and simply normalizes such that the max sharpness value is 1 and the min
% is 0.
%
%IN
% global_sharpness: (1 x n_image) array with the sharpness metric for each
% entire image
%
% window_sharpness: (n_window x n_image) array with the sharpness metric
% for each window specified
%
% OUT 
% global_normalized: (1 x n_image) array with the sharpness metric for each
% entire image normalized
%
% window_sharpness: (n_window x n_image) array with the sharpness metric
% for each window normalized with respect to that own window's sharpness
%
%
% R. A. Manzk 12/08/2020
%% begin the function
    min_global = min(global_sharpness,[],2);
    max_global = max(global_sharpness,[],2);
    global_normalized = (global_sharpness - min_global)./(max_global-min_global);

    min_window = min(window_sharpness,[],2);
    max_window = max(window_sharpness,[],2);
    window_normalized = (window_sharpness - min_window)./(max_window-min_window);
    end