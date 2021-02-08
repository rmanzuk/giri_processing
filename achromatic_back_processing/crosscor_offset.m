function calculated_offset = crosscor_offset(shifted_im,fixed_im,edges)
% this function calcualtes the offset in pixels between two images using
% their normalized cross correlation
%
% IN 
%
% shifted_im: m x n array representing single_channel input image as a
% double. This is the one that is considered to be shifted relative to the
% other
%
% fixed_im: m x n array representing single_channel input image as a
% double. This is the one that is considered to be fixed or unmoving.
%
% edges: binary flag indicating if edges of the image should be used as
% opposed to the raw image. 1 for edge, 0 for plain image.
%
% OUT
%
% calculated_offset: 2x1 array with the x and y offsets between the images
% to be applied as a shift in pixel space for alignment
    %% begin the function
    % if we just want edges, find them
    if edges
        shifted_im = edge(shifted_im, 'canny');
        fixed_im = edge(fixed_im, 'canny');
    else
        % do nothing
    end
    % if the image is flat, we can't do anything
    if std(fixed_im(:)) == 0 || std(shifted_im(:)) == 0
        calculated_offset = [NaN,NaN];
    else
        
        % calculate cross corr
        nxc = normxcorr2(shifted_im,fixed_im);

        % index the peak
        [~,index_max] = max(abs(nxc(:)));
        [ypeak,xpeak] = ind2sub(size(nxc),index_max(1));

        % back into image coordinates
        calculated_offset = [(ypeak-size(shifted_im,1)),(xpeak-size(shifted_im,2))];
    end
    
end