function [superpix_training_data] = get_superpix_training(input_ims,label_mats,labeled_pix)

    %% begin the function
    % empty matrix for the training data to end up in
    superpix_training_data = [];
    
    % also good to know the fields in the structs
    fields = fieldnames(label_mats);

    % loop through the number of classes 
    for i = 1:numel(labeled_pix.(fields{1}))
        % and loop through the number of input images
        for j = 1:numel(fields)
            % grab this selection of superpixels
            this_label_mat = label_mats.(fields{j});
            % grab the pixels labeled for this class
            these_pix = labeled_pix.(fields{j}){i};
            % get the index cell array for the superpixels
            idx = label2idx(this_label_mat);
            % and just select this image
            this_im = input_ims.(fields{j});
            % empty array to take these training data
            these_training_data = zeros(numel(these_pix),size(this_im,3)*2);
            % now we gotta loop through the super pixels for each class
            for k = 1:numel(these_pix)
                % unfortunately, I think we need to loop through each channel.
                % Maybe look for a faster solution in the future
                for j = 1:size(this_im,3)
                    % grab linear indeces based upon the channel number
                    inds = idx{these_pix(k)} +(j-1)*size(this_im,1)*size(this_im,2);
                    % use inds to extract values
                    pixel_values = this_im(inds);
                    % put the mean and standard deviation into place in the stats
                    % mat
                    these_training_data(k,j) = mean(pixel_values);
                    these_training_data(i,j+size(this_im,3)) = std(pixel_values);
                end
            end
            % add in a final column for class
            these_training_data(:,end+1) = i;
            superpix_training_data = [superpix_training_data;these_training_data];
        end

    end
end