function [sharpest_image, sharpest_loc] = select_sharpest(input_dir, global_sharpness, window_sharpness, method, n_zsteps, z_step_size)
% This function takes the sharpness metric array outputs of the
% stack_sharpness function and outputs the sharpest image. The user can
% specify whether that image is selected based upon the whole image or
% window sharpness or both.
%
%IN
% input_dir: pathname string for the folder containing the stack of images
% for a channel
%
% global_sharpness: (1 x n_image) array with the sharpness metric for each
% entire image
%
% window_sharpness: (n_window x n_image) array with the sharpness metric
% for each window specified
% 
% method: 'whole' to use the whole image sharpness, 'windows' to use the
% mean window sharpness, 'both' to use both additively
%
% n_zsteps: number of z_steps above and below the central (0) image plane.
% For example, if you have 31 total images with a central image, 15
% above, and 15 below, this input would be 15.
% 
% z_step_size: size of each z_step in microns.
%
% OUT 
% sharpest_image: uint8 or 16 (depending on input bit depth) 3D array with
% the most focused image for that chanel.
%
% sharpest_loc: location of the sharpest image in microns relative to the
% central (0) image plane
%
%
% R. A. Manzk 11/20/2020
    %% begin the function
    % set up the directory and file reading
    file_pattern = fullfile(input_dir, '*.tif');
    tifs = dir(file_pattern);
    base_names = natsortfiles({tifs.name});
    
    % depending on method, take the max sharpness value
    if strcmp(method,'whole')
        [~,sharpest_ind] = max(global_sharpness);
    elseif strcmp(method,'windows')
        [~,sharpest_ind] = max(mean(window_sharpness));
    else
        [~,sharpest_ind] = max(mean(window_sharpness)+global_sharpness);
    end
    
    % use the sharpness values as an index to read the stack and grab
    % sharpest image
    sharpest_name = fullfile(input_dir, base_names{sharpest_ind});
    sharpest_image = imread(sharpest_name);
    
    % and based on that index, where are we geometrically
    sharpest_loc = (sharpest_ind - n_zsteps) * z_step_size;
end