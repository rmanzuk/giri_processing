function [global_sharpness,window_sharpness] = stack_sharpness(input_dir,n_windows,window_size,n_zsteps,z_step_size)
% THis function goes through a single color channel's z-stack (organized in
% a folder and takes an overall sharpness metric for each as well as
% sharpness for a number of smaller windows specified
% 
% IN
% input_dir: pathname string for the folder containing the stack of images
% for a channel
%
% n_windows: the number of smaller windows of which you would like to judge
% sharpness
%
% window_size: number of pixels for the window size. Should be a single
% number as windows will be square
% 
% n_zsteps: number of z_steps above and below the central (0) image plane.
% For example, if you have 31 total images with a central image, 15
% above, and 15 below, this input would be 15.
% 
% z_step_size: size of each z_step in microns.
%
% OUT
% global_sharpness: (1 x n_image) array with the sharpness metric for each
% entire image
%
% window_sharpness: (n_window x n_image) array with the sharpness metric
% for each window specified
%
% R. A. Manzuk 11/19/2020
    %% begin the function
    % set up the directory and file reading
    file_pattern = fullfile(input_dir, '*.tif');
    tifs = dir(file_pattern);
    base_names = natsortfiles({tifs.name});
    n_images = numel(base_names);

    % just load the frist picture to take some random windows for sharpness checks
    sample_im_name = fullfile(input_dir, base_names{1});
    sample_im = imread(sample_im_name);
    window_coords = round(rand(n_windows,2).*size(sample_im));
%     imshow(sample_im)
%     fprintf('Plese select 10 points for sharpness windows %u\n', n_windows)
%     window_coords = ginput(10);
%     window_coords = round(window_coords);
    

    % simple array counting number of images for plotting
    im_count = [0:(numel(tifs)-1)];
    
    % set up arrays to catch sharpness metrics
    global_sharpness = zeros(1,numel(tifs));
    window_sharpness = zeros(n_windows, numel(tifs));
    %iterate through images
    for i = 1:numel(tifs)
        tic
        % read in image
        fprintf('Now reading image %u of %u\n', [i,numel(tifs)]);
        full_file_name = fullfile(input_dir, base_names{i});
        this_im = imread(full_file_name);
        this_im = im2double(this_im);
        % estimate of sharpness as the sum of all gradient norms/number of
        % pixels
        [grad_x,grad_y] = gradient(this_im);
        grad_norms = sqrt(grad_x.*grad_x+grad_y.*grad_y);
        global_sharpness(i) = sum(sum(grad_norms))./(numel(grad_norms));
        % then go and get the sharpenss just in the windows with same metric
        for j = 1:n_windows
            row_ind1 = window_coords(j,1);
            col_ind1 = window_coords(j,2);
            row_ind2 = row_ind1 + window_size;
            col_ind2 = col_ind1 + window_size;

            if row_ind2 > size(this_im,1)
                row_ind2 = size(this_im,1);
            end

            if col_ind2 > size(this_im,2)
                col_ind2 = size(this_im,2);
            end
            window_grad = grad_norms(row_ind1:row_ind2,col_ind1:col_ind2);
            window_sharpness(j,i) = sum(sum(window_grad))./(numel(window_grad));
        end

        % and plot how it's going
        subplot(2,1,1)
        plot((im_count(1:i)-n_zsteps).*z_step_size,global_sharpness(1:i),'DisplayName','global sharpness')
        ylabel('Sharpness')
        title('Whole image')
        subplot(2,1,2)
        plot((im_count(1:i)-n_zsteps).*z_step_size,window_sharpness(:,1:i),'LineStyle',':')
        xlabel('Z location [um]')
        ylabel('Sharpness')
        title('Windows')
        drawnow

        toc


    end
    end