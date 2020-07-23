function [] = align_stack(fixed, moving, rest_of_stack, ext, output_directory)
% This function uses feature detection and matching to re-align GIRI
% image stacks when they get off-kilter. Before using, split the stack into
% separate folders at the point of misalignment, and load into matlab the
% images just above and below that point. Then just press go and follow the
% prompts
%
% IN
% fixed: 3D (3 channel) image of the slice just above the point of misalignment. This
% is the reference point to which the lower images in the stack will be
% aligned
%
% moving: 3D (3 channel) image of the slice just below the point of
% misalignment. We'll align this one to the 'fixed' slice and use that
% transformation for the rest of the lower stack
%
% rest_of_stack: string with the file path for the folder containing all of
% the images below the shift point to be aligned.
% e.g. '/Users/rmanzuk/Desktop/code/giri_processing/stack_alignment/lower'
% 
% ext: string for the extension for the image files in the rest_of_stack
% folder with a star in front.... e.g. '*.tif'
%
% output_directory: string with the file path for the folder where you
% would like the aligned images to go. 
% e.g. '/Users/rmanzuk/Desktop/code/giri_processing/stack_alignment/lower_aligned'
%
%
% OUT 
% nothing.....it just puts the aligned images where they belong
% R. A. Manzuk, 07/23/2020
    % begin the function
    % identify SURF features from the two images to be aligned
    fixed_points_SURF = detectSURFFeatures(rgb2gray(fixed));
    moving_points_SURF = detectSURFFeatures(rgb2gray(moving));

    % extract SURF descriptors of the features
    [fixed_features_SURF, fixed_valid_SURF] = extractFeatures(rgb2gray(fixed),fixed_points_SURF);
    [moving_features_SURF, moving_valid_SURF] = extractFeatures(rgb2gray(moving),moving_points_SURF);

    % and match up the SURF features at many thresholds
    thresholds_surf = [0.1,0.2,0.3,0.4];
    feature_pairs_SURF = {};
    fixed_matched_SURF = {};
    moving_matched ={};
    for i = 1:length(thresholds_surf)
        feature_pairs_SURF{i} = matchFeatures(fixed_features_SURF,moving_features_SURF,'MatchThreshold',thresholds_surf(i));
        fixed_matched_SURF{i} = fixed_valid_SURF(feature_pairs_SURF{i}(:,1),:);
        moving_matched_SURF{i} = moving_valid_SURF(feature_pairs_SURF{i}(:,2),:);
    end

    % identify Harris features from the two images to be aligned
    fixed_points_harris = detectHarrisFeatures(rgb2gray(fixed));
    moving_points_harris = detectHarrisFeatures(rgb2gray(moving));

    % extract Harris descriptors of the features
    [fixed_features_harris, fixed_valid_harris] = extractFeatures(rgb2gray(fixed),fixed_points_harris);
    [moving_features_harris, moving_valid_harris] = extractFeatures(rgb2gray(moving),moving_points_harris);

    % and match up the Harris features at many thresholds
    thresholds_harris = [1,5,10,20];
    feature_pairs_harris = {};
    fixed_matched_harris = {};
    moving_matched_harris ={};
    for i = 1:length(thresholds_surf)
        feature_pairs_harris{i} = matchFeatures(fixed_features_harris,moving_features_harris,'MatchThreshold',thresholds_harris(i),'MaxRatio',0.8);
        fixed_matched_harris{i} = fixed_valid_harris(feature_pairs_harris{i}(:,1),:);
        moving_matched_harris{i} = moving_valid_harris(feature_pairs_harris{i}(:,2),:);
    end


    % show the SURF matches
    surf_fig = figure(1);
    subplot(2,2,1)
    showMatchedFeatures(fixed,moving,fixed_matched_SURF{1},moving_matched_SURF{1});
    str1 = sprintf('threshold = %.1f',thresholds_surf(1));
    title(str1)
    subplot(2,2,2)
    showMatchedFeatures(fixed,moving,fixed_matched_SURF{2},moving_matched_SURF{2});
    str2 = sprintf('threshold = %.1f',thresholds_surf(2));
    title(str2)
    subplot(2,2,3)
    showMatchedFeatures(fixed,moving,fixed_matched_SURF{3},moving_matched_SURF{3});
    str3 = sprintf('threshold = %.1f',thresholds_surf(3));
    title(str3)
    subplot(2,2,4)
    showMatchedFeatures(fixed,moving,fixed_matched_SURF{4},moving_matched_SURF{4});
    str4 = sprintf('threshold = %.1f',thresholds_surf(4));
    title(str4)

    sgtitle('SURF feature matches')

    % show the Harris matches
    harris_fig = figure(2);
    subplot(2,2,1)
    showMatchedFeatures(fixed,moving,fixed_matched_harris{1},moving_matched_harris{1});
    str1 = sprintf('threshold = %.1f',thresholds_harris(1));
    title(str1)
    subplot(2,2,2)
    showMatchedFeatures(fixed,moving,fixed_matched_harris{2},moving_matched_harris{2});
    str2 = sprintf('threshold = %.1f',thresholds_harris(2));
    title(str2)
    subplot(2,2,3)
    showMatchedFeatures(fixed,moving,fixed_matched_harris{3},moving_matched_harris{3});
    str3 = sprintf('threshold = %.1f',thresholds_harris(3));
    title(str3)
    subplot(2,2,4)
    showMatchedFeatures(fixed,moving,fixed_matched_harris{4},moving_matched_harris{4});
    str4 = sprintf('threshold = %.1f',thresholds_harris(4));
    title(str4)

    sgtitle('Harris feature matches')

    % maybe the user isn't happy with any of those options, so ask if they'd
    % like to try and see anything else
    unhappy = 1;
    while unhappy
        prompt = 'Would you like to see any method with a different threshold? (Y/N) \n';
        more_match = input(prompt,'s');
        if strcmpi(more_match, 'N')
            unhappy = 0;
        elseif strcmpi(more_match, 'Y')
            prompt2 = 'Which feature detection method? (1 for SURF, 2 for Harris) \n';
            method = input(prompt2);
            prompt3 = 'Which matching threshold would you like to use? (input number) \n';
            thresh = input(prompt3);

            if method == 1
                feature_pairs = matchFeatures(fixed_features_SURF,moving_features_SURF,'MatchThreshold',thresh);
                fixed_matched = fixed_valid_SURF(feature_pairs(:,1),:);
                moving_matched = moving_valid_SURF(feature_pairs(:,2),:);
                figure() 
                showMatchedFeatures(fixed,moving,fixed_matched,moving_matched);
                string = sprintf('SURF threshold = %.1f',thresh);
                title(string)
            elseif method == 2
                feature_pairs = matchFeatures(fixed_features_harris,moving_features_harris,'MatchThreshold',thresh,'MaxRatio',0.8);
                fixed_matched = fixed_valid_harris(feature_pairs(:,1),:);
                moving_matched = moving_valid_harris(feature_pairs(:,2),:);
                figure()
                showMatchedFeatures(fixed,moving,fixed_matched,moving_matched);
                string = sprintf('Harris threshold = %.1f',thresh);
                title(string)
            else
                disp('Try again')
            end   
        else
            disp('Tell Ryan you are not happy')
        end
    end

    close all

    % now that the user has settled on good feature matching, ask them for
    % parameters, and execute the transformation
    prompt4 = 'For final alignment, which method? (1 for SURF, 2 for Harris) \n';
    method2 = input(prompt4);
    prompt5 = 'Which matching threshold would you like to use for final aligment? (input number) \n';
    thresh2 = input(prompt5);

    % execute the matching with user-desired parameters
    if method2 == 1
        feature_pairs = matchFeatures(fixed_features_SURF,moving_features_SURF,'MatchThreshold',thresh2);
        fixed_matched = fixed_valid_SURF(feature_pairs(:,1),:);
        moving_matched = moving_valid_SURF(feature_pairs(:,2),:);
    elseif method2 == 2
        feature_pairs = matchFeatures(fixed_features_harris,moving_features_harris,'MatchThreshold',thresh2,'MaxRatio',0.8);
        fixed_matched = fixed_valid_harris(feature_pairs(:,1),:);
        moving_matched = moving_valid_harris(feature_pairs(:,2),:);
        figure()
    else
    end


    % solve for the image transform that gives proper alignment
    transform = estimateGeometricTransform(moving_matched,fixed_matched, 'similarity');
    % transform the moving image given that transform

    test_transformed = imwarp(moving,transform,'OutputView',imref2d(size(fixed)));

    % let's see how we did, display the entire images
    figure(1)
    imshowpair(fixed,test_transformed)
    title('Whole alignment')

    % and look at some smaller areas
    % select random coordinates for zoomed-in windows to really test alignment
    window_coords = (round(rand(4,2).* [size(fixed,2)*(6/8),size(fixed,1)*(6/8)])) + [size(fixed,2)/8,size(fixed,1)/8];

    % extract the neighborhoods
    fixed_neighborhood = get_pixel_neighborhood(fixed,window_coords,100);
    transformed_neighborhood = get_pixel_neighborhood(test_transformed,window_coords,100);

    % show them
    figure(2)
    subplot(2,4,1)
    imshow(fixed_neighborhood(:,:,:,1))
    title('1 - upper slice')
    subplot(2,4,2)
    imshow(fixed_neighborhood(:,:,:,2))
    title('2 - upper slice')
    subplot(2,4,3)
    imshow(fixed_neighborhood(:,:,:,3))
    title('3 - upper slice')
    subplot(2,4,4)
    imshow(fixed_neighborhood(:,:,:,4))
    title('4 - upper slice')
    subplot(2,4,5)
    imshow(transformed_neighborhood(:,:,:,1))
    title('1 - lower slice')
    subplot(2,4,6)
    imshow(transformed_neighborhood(:,:,:,2))
    title('2 - lower slice')
    subplot(2,4,7)
    imshow(transformed_neighborhood(:,:,:,3))
    title('3 - lower slice')
    subplot(2,4,8)
    imshow(transformed_neighborhood(:,:,:,4))
    title('4 - lower slice')

    sgtitle('Zoomed in tiles')
    
    % ask the user if they are happy
    prompt = 'Are you satisfied with the alignment? (Y/N) \n';
    satisfaction = input(prompt,'s');

    if strcmpi(satisfaction, 'Y')
        disp('Well that is terrific.....processing rest of stack now');
        file_pattern = fullfile(rest_of_stack, ext);
        input_files = dir(file_pattern);
        ordered_files = natsortfiles({input_files.name});
        for i = 1:numel(ordered_files)
            file_to_read = fullfile(rest_of_stack, ordered_files{i});
            file_to_write = fullfile(output_directory, ordered_files{i});
            [~,image_name,~] = fileparts(file_to_read);
            fprintf('Now reading %s\n', file_to_read);
            image = imread(file_to_read);
            transformed_image = imwarp(image,transform,'OutputView',imref2d(size(fixed)));
            imwrite(transformed_image,file_to_write,'tiff')
        end
        return
    elseif strcmpi(satisfaction, 'N')
        % maybe cycle through the alignment with another method?
    else
    end
end