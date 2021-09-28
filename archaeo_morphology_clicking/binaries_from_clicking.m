function [] = binaries_from_clicking(inner_data,outer_data,output_dir)
% This function takes the manually clicked data from an archaeo clicking
% stack and outputs binary tiffs for the resulting model.
%
% IN
% outer_data: 1xn_branches cell array containing the cell arrays for outer circles 
% created during data collection. could be from clicking or automated
% tracing.
%
% inner_data: 1xn_branches cell array containing the cell arrays for inner circles 
% created during data collection. could be from clicking or automated
% tracing.
% 
% output_dir: string with the path to a folder where you would like the
% binary tiffs to be saved.
%
% OUT
%
%
% R. A. Manzuk 09/28/2021
    %% begin the function
    % find the number of slices in total
    n_slices = max(cellfun(@numel,outer_data));

    % make sure all individual branch cell arrays are same size...just in
    % case
    for i = 1:numel(outer_data)
        if numel(outer_data{i}) < n_slices
            add_empty = cell(1,n_slices - numel(outer_data{i}));
            outer_data{i}(numel(outer_data{i})+1:n_slices) = add_empty;
        end
        if numel(inner_data{i}) < n_slices
            add_empty = cell(1,n_slices - numel(inner_data{i}));
            inner_data{i}(numel(inner_data{i})+1:n_slices) = add_empty;
        end
    end

    % how big is the model in the x-y plane?
    % set up empty array, and just go through the cell array and stack all 
    % points into one 2-column matrix
    all_xy = [];
    for i = 1:numel(outer_data)
        for j = 1:n_slices
            if ~isempty(outer_data{i}{j})
                all_xy = [all_xy; outer_data{i}{j}(:,1:2)];
            end
        end
    end

    % then just gather the max, min of those xy coord and use to define size of
    % image plane
    min_xy = min(all_xy);
    max_xy = max(all_xy);

    im_dims = ceil(max_xy-min_xy);

    % now just go through each slice and make the image
    for i = 1:n_slices
       this_im = zeros(flip(im_dims)); 
       for j = 1:numel(outer_data)
          % make a mask of outer data if exists
          if ~isempty(outer_data{j}{i})
              branch_mask = poly2mask(outer_data{j}{i}(:,1)-min_xy(1),outer_data{j}{i}(:,2)-min_xy(2),im_dims(2),im_dims(1));
          end
          % make a mask of inner data if exists
          if ~isempty(inner_data{j}{i})
              osc_mask = poly2mask(inner_data{j}{i}(:,1)-min_xy(1),inner_data{j}{i}(:,2)-min_xy(2),im_dims(2),im_dims(1));
          end
       % apply the branch mask
       this_im(branch_mask == 1) = 1;
       % apply the osculum mask
       this_im(osc_mask == 1) = 0;
       end
       % write the image
       im_name = sprintf('binary_%d.tiff',i);
       file_name = fullfile(output_dir, im_name);
       imwrite(this_im,file_name);
    end
end