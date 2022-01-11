function [multichannel_im] = concat_norm_im(varargin)
% This function takes in any number of channels as doubles, maximally
% stretches their histograms and concatenates them
%
% IN: 
% varargin: inputs should be each individual image channel (2D double
% matrix) separated by commas, in the order for concatenation.
%
% OUT:
% multichannel_im: 3D double matrix containing all channels, stretched and
% concatenated
%
% Ryan A. Manzuk 11/15/2021
    %%
    % number of channels will be defined by number or inputs
    num_channels = nargin;
    % set up empty array for the final image
    multichannel_im = zeros(size(varargin{1},1),size(varargin{1},2),nargin);

    % loop through, normalize, and place each channel
    for i = 1:nargin
        multichannel_im(:,:,i) = (varargin{i} - min(varargin{i}(:))) * (1./max(varargin{i}(:)));
    end
end
