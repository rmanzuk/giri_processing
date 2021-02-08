function [im_MSE] = best_blur(sharp_im,blurred_im,sigma,net)
% This function simultes the blur between two images with gaussian filters 
% of various sizes (sizes based on input sigma) after quantifying image 
% noise with a CNN. It outputs the MSE between the simulated blur and
% actual blurred image for each sigma. 
%
% IN: 
% sharp_im: mxn sharp image, as a grayscale double. This image will undergo
% the simulated blurring to try to match the blurred_im.
%
% blurred_im: mxn blurred image, as a grayscale double.
% 
% sigma: 1D vector with the sigma values to be used for the gaussians to
% simulate blur. Somthing like [0.1:0.1:10] would probably work well.
%
% net: denoising neural net input as a SeriesNetwork object. Try 
% net = denoisingNetwork('DnCNN');
%
% OUT:
% im_MSE: 1D vector containing the mean sqared error between each simulated
% blurred image and the actual blurred image for every sigma gaussian used.
% The min of this vector would indicate the best blur simulation.
%
% Ryan A. Manzuk 01/28/2021
%%
    denoised_sharp = denoiseImage(sharp_im,net);
    noise = sharp_im - denoised_sharp;
    noise_var = var(noise(:));

    im_MSE = zeros(1,length(sigma));
    for i = 1:length(sigma)
        simulated_blur = imgaussfilt(denoised_sharp,sigma(i));
        simulated_blur_noisy = imnoise(simulated_blur,'gaussian',0,noise_var);
        im_MSE(i) = (sum((simulated_blur_noisy - blurred_im).^2,'all'))/numel(blurred_im);
    end
end