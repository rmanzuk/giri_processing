function [pc_image, pc_loadings, eigens, covariance_mat] = imPCA(input_im)
% This function takes an image and does the reshaping and normalization
% necessary to do a PCA on its color data and returns an image of the same
% dimensions where pixel values are now pc scores, and the channels are the
% different principal components
%
% IN
% input_im: multi-channel image in any format (double, 16-bit, 8-bit, etc.)
%
% OUT
%
% pc_image: image of the same dimensions as the input image where each
% channel is a principal component and the pixel values are the pc scores.
% The values are such that the maximum score is 1 and the minimum is 0, but
% there is no normalization between the principal components. 
%
% pc_loadings: n_channel x n_channel square matrix for the loadings of each
% of the image channels on the principal components.
%
% eigens: eigenvalues of the principal components. Can be used to derive
% the pct of variance explained by each or for later
% weighting/normalization
%
% covariance_mat: diagonal covariance matrix for all of the image channels
%
% R. A. Manzuk 03/26/2021
    %% begin the function
    % need to know those sizes
    [dim1, dim2, n_channels] = size(input_im);

    % reshape image into n_pixel x n_channel 2d array
    reshaped_im = reshape(input_im, dim1 * dim2, n_channels);

    % normalize such that all channels are 0-center win standard deviation
    % of 1
    zero_center = reshaped_im - mean(reshaped_im);
    normalized_im = zero_center./std(zero_center);
    
    % might be nice to know the covariance
    covariance_mat = cov(normalized_im);
    
    % run the pca
    [pc_loadings, pc_scores, eigens] = pca(normalized_im);

    % bring the scores back between 0 and 1
    pc_scores = pc_scores - min(pc_scores);
    pc_scores = pc_scores./max(pc_scores);

    % reshape back to image
    pc_image = reshape(pc_scores, dim1, dim2, n_channels);
end