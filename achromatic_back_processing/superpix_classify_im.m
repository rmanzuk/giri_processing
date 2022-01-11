function [classified_im] = superpix_classify_im(model,to_classify,label_mat)
% this function takes a trained model and uses it to classify an entire
% image after superpixel oversegmentation.
%
% IN:
%
% model: trained model (I've used svms and random forests) in the class
% given by matlab.
%
% to_classify: Double image with the same number of channels as the those
% used for training of the model.
%
% label_mat: Superpixel label matrix, given by the superpixel algorithm in
% matlab.
%
% OUT: 
%
% classified_im: single channel image with the same row and collumn
% dimensions as the input image where all pixels are given the index
% correstponding to their assigned class.
%
% R. A. Manzuk 11/30/2021
    %% begin the function
    % good to know the number of channels we have
    n_channels = size(to_classify,3);

    % get the index cell array for the superpixels
    idx = label2idx(label_mat);

    % empty matrix for the superpixel stats
    superpix_stats = zeros(max(label_mat,[],'all'),n_channels*2);

    % loop through each superpixel and grab the stats
    for i = 1:max(label_mat,[],'all')
        % unfortunately, I think we need to loop through each channel.
        % Maybe look for a faster solution in the future
        for j = 1:n_channels
            % grab linear indeces based upon the channel number
            inds = idx{i} +(j-1)*size(to_classify,1)*size(to_classify,2);
            % use inds to extract values
            pixel_values = to_classify(inds);
            % put the mean and standard deviation into place in the stats
            % mat
            superpix_stats(i,j) = mean(pixel_values);
            superpix_stats(i,j+n_channels) = std(pixel_values);
        end
    end

    % use those stats to classify each pixel
    predicted_supers = predict(model,superpix_stats);

    if iscell(predicted_supers)
        predicted_supers = cell2mat(predicted_supers);
        predicted_supers = str2num(predicted_supers);
    end

    % now we have to put those predictions into a classified im
    classified_im = zeros(size(to_classify,1),size(to_classify,2));
    % loop through each superpixel
    for i = 1:numel(idx)
        % put the prediction in the indeces
        classified_im(idx{i}) = predicted_supers(i);
    end

end