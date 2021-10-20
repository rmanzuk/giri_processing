function [positive_glcm_stats, negative_glcm_stats] = masked_glcm(input_im, mask)
% this function takes an input image with a binary, logical mask that
% separates two classes and gives the glcm stats for each class (positive
% ids and negative ids)
%
% IN
% input_im: 2D (single channel) double image
%
% mask: logical mask of same dimensions as input_im where true represents a
% positive id for the desired class an false represents an negative id
%
% OUT
%
% positive_glcm_stats: structure containing the 4 GLCM statistics for
% positively identified regions
%
% negative_glcm_stats:structure containing the 4 GLCM statistics for
% negatively identified regions
%
% R. A. Manzuk 03/26/2021
    %% begin the function
    
    % set the negative id pixels to Nan
    positive_im = input_im;
    positive_im(~mask) = NaN;
    
    % do the GLCM stuff
    positive_glcm = graycomatrix(positive_im);
    positive_glcm_stats = graycoprops(positive_glcm);
    
    % set positive id pixels to nan
    negative_im = input_im;
    negative_im(mask) = NaN;
    
    % do the GLCM stuff
    negative_glcm = graycomatrix(negative_im);
    negative_glcm_stats = graycoprops(negative_glcm);

end