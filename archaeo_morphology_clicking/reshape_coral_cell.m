function [reshaped_outers] = reshape_coral_cell(initial_outers)
% this function takes the 1xn_slices cell array from the automated tracing
% output and reshapes it to 1xn_branches array that can be used in the
% measurement pipeline.

% IN
% initial_outers: 1xn_slice cell array containing the edge tracings from
% the automated tracing pipeline, could be cleaned up with spurious edges
% removed.
%
% OUT
% reshaped_outers: 1xn_branch cell array containing the edge tracings from
% the automated tracing pipeline, reshaped for input to measurement
% functions.
%
% R. A. Manzuk 12/21/2020
    %% begin the function
    % get some initial info
    total_edges = max(cellfun(@numel,initial_outers));
    
    
    %figure out which edges are not present in which slices
    not_empty_ind = zeros(size(initial_outers,2),total_edges);
    for i = 1:numel(initial_outers)
        not_empty_ind(i,1:numel(initial_outers{i})) = ~cellfun(@isempty,initial_outers{i});
    end

    % set up the reshaped array
    reshaped_outers = cell(1,total_edges);
    
    % fill in the reshaped array with all of the info
    for i = 1:numel(reshaped_outers)
        reshaped_outers{i} = cell(1,size(initial_outers,2));
        for k = 1:numel(reshaped_outers{i})
            % only look to fill spot if there are data
            if not_empty_ind(k,i)
                reshaped_outers{i}{k} = initial_outers{k}{i};
            else
            end
        end
    end
end