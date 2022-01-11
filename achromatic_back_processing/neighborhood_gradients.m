function [Gmag_neigh,Gdir_neigh] = neighborhood_gradients(input_im,rad)
% This function takes an input image and assesses its gradient magnitude
% and direction. It then gives back images that are the local averages at
% every pixel for the gradients within a circular neighborhood of specified
% radius. If the input image is multichannel, the function reduces it to
% its first PC values.
%
% IN
% input_im: 2d or 3d matrix with double image (can be multichannel or
% grayscale)
%
% rad: radius of the circular neighboorhood you would like to average
% within for the gradient magnitude and direction images
%
% OUT
%
% Gmag_neigh: 2d matrix with the reasulting averaged local gradient
% magnitudes within the surrounding neighborhood.
%
% Gdir_neigh: 2d matrix with the reasulting averaged local gradient
% directions within the surrounding neighborhood.
%
% R. A. Manzuk 11/23/2021
    %% begin the function
    % if we have more than 1 channel, we need to reduce, so let's just use a
    % pca
    if size(input_im,3) > 1
        % do the PCA with a prewritten function
        pc_im = imPCA(input_im);
        % just take the first channel/PC for max variance
        to_process = pc_im(:,:,1);
    else
        % if we just have 1 channel, we don't need to do anything
        to_process = input_im;
    end
    
    % then just extract the gradient magnitude and direction
    [Gmag,Gdir] = imgradient(to_process);
    
    % make a circular filter that has the radius specified
    filt = fspecial('disk',rad);
    
    % filter the gradient magnitudes and directions with that filter to make
    % images representing the local averages of those things
    
    Gmag_neigh = conv2(Gmag,filt,'same');
    Gdir_neigh = conv2(Gdir,filt,'same');
end
