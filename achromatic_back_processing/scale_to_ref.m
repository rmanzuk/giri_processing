function [scaled_im,transform] = scale_to_ref(reference_channel_im, im_to_scale)
% This function takes one channel image as the reference image and uses
% feature detection and matching to scale the other 
%
%IN
%
%
% OUT 
%
%
% R. A. Manzuk 11/20/2020
    %% begin the function
    reference_features = detectSURFFeatures(im2double(reference_channel_im));
    to_scale_features = detectSURFFeatures(im2double(im_to_scale));
    
    % extract SURF descriptors of the features
    [reference_features, reference_valid] = extractFeatures(im2double(reference_channel_im),reference_features);
    [to_scale_features, to_scale_valid] = extractFeatures(im2double(im_to_scale),to_scale_features);
    % and match up the SURF features at many thresholds
    thresholds = [0.1,0.2,0.3,0.4];
    feature_pairs = {};
    reference_matched = {};
    to_scale_matched ={};
    for i = 1:length(thresholds)
        feature_pairs{i} = matchFeatures(reference_features,to_scale_features,'MatchThreshold',thresholds(i));
        reference_matched{i} = reference_valid(feature_pairs{i}(:,1),:);
        to_scale_matched{i} = to_scale_valid(feature_pairs{i}(:,2),:);
    end
    
    % show the SURF matches
    surf_fig = figure(1);
    subplot(2,2,1)
    showMatchedFeatures(reference_channel_im, im_to_scale,reference_matched{1},to_scale_matched{1});
    str1 = sprintf('threshold = %.1f',thresholds(1));
    title(str1)
    subplot(2,2,2)
    showMatchedFeatures(reference_channel_im, im_to_scale,reference_matched{2},to_scale_matched{2});
    str2 = sprintf('threshold = %.1f',thresholds(2));
    title(str2)
    subplot(2,2,3)
    showMatchedFeatures(reference_channel_im, im_to_scale,reference_matched{3},to_scale_matched{3});
    str3 = sprintf('threshold = %.1f',thresholds(3));
    title(str3)
    subplot(2,2,4)
    showMatchedFeatures(reference_channel_im, im_to_scale,reference_matched{4},to_scale_matched{4});
    str4 = sprintf('threshold = %.1f',thresholds(4));
    title(str4)
    
    sgtitle('SURF feature matches')
    
    % ask the user which threshold is best
    prompt = 'Which matching threshold would you like to use? (input number) \n';
    thresh = input(prompt);
    
    % execute the matching with user-desired parameters
    feature_pairs = matchFeatures(reference_features,to_scale_features,'MatchThreshold',thresh);
    reference_matched = reference_valid(feature_pairs(:,1),:);
    to_scale_matched = to_scale_valid(feature_pairs(:,2),:);
    
    % solve for the image transform that gives proper alignment
    transform = estimateGeometricTransform(reference_matched,to_scale_matched, 'similarity');
    
    % transform the moving image given that transform
    scaled_im = imwarp(reference_channel_im, im_to_scale,'OutputView',imref2d(size(reference_channel_im)));
end