Efunction [] = transverse_resample_stack(input_folder, ext, output_folder, output_naming, scale_ratio, resampling_dim)
% This function takes and existing grinder image stack and resamples it
% along one of the transverse axes to give new, edge-on photos
%
% IN
% input_folder: full pathway to the folder containing the images for the
% original grinder image stack
%
% ext: file extenstion for the image type in *. form. Example '*.tiff'
%
% output_folder: full pathway to the folder where you would like the new
% images to go.
%
% output_naming: string with the image name you would like to be at the
% front of each new image file. Example: for images named RM_transverse01,
% RM_transvers02, etc...this variable would be 'RM_transverse'
%
% scale_ratio: ratio between the vertical slice spacing and the per_pixel
% resolution for the stack. Example: if your images have a resolution of
% 5.4 um/pixel and your stack has 20 um spacing, this variable would be
% 20/5.4
%
% resampling_dim: dimension along which you would like to resample images
% to give transverse section. 1 to sample row-wise, 2 to sample
% column-wise.
%
% OUT
%
% R. A. Manzuk 10/08/2020
%% begin the function
    % what is the non-resampling dimension?
    if resampling_dim == 1
        other_dim = 2;
    else
        other_dim = 1;
    end

    % set up the input folder file pattern
    file_pattern = fullfile(input_folder, ext);
    tifs = dir(file_pattern);
    base_names = natsortfiles({tifs.name});
    n_images = numel(base_names);

    % just load the frist picture to set up dimensions
    sample_im_name = fullfile(input_folder, base_names{1});
    sample_im = imread(sample_im_name);


    for i = 1:1%size(sample_im,resampling_dim)
        tic
        %set up to make the new image
        fprintf('Now making image %u of %u\n', [i,size(sample_im,resampling_dim)]);
        new_im = zeros(round(n_images*scale_ratio),size(sample_im,other_dim),3);
        count = 0;
        for j = 1:numel(base_names)
            % going to need to grab every image 
            full_file_name = fullfile(input_folder, base_names{j});
            this_im = imread(full_file_name);
            %bookkeeping 
            depth_difference = round(j*scale_ratio) - count;
            % extract the proper line to be placed on the current image
            if resampling_dim == 1
                this_line = this_im(i,:,:);
                to_place = repmat(this_line,depth_difference,1,1);
                new_im((count+1):(count+depth_difference),:,:) = to_place;
            else
                this_line = this_im(:,i,:);
                to_place = repmat(this_line,1,depth_difference,1);
                new_im(:,(count+1):(count+depth_difference),:) = to_place;
            end
            count = count + depth_difference;
        end
        output_file_name = fullfile(output_folder, strcat(output_naming,string(i)));
        imwrite(new_im,output_file_name);
        toc
    end
end

