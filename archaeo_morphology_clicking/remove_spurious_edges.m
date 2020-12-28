function [cleaned_outers] = remove_spurious_edges(initial_outers,nslice_cutoff)
% This function goes through and removes edges that don't appear on enough
% slices, and so must be considered spurious from the automated tracing
% output.
%
% IN
% initial_outers: 1xn_slice cell array containing the edge tracings from
% the automated tracing pipeline
%
% nslice_cutoff: Number of slices appearing on, below which an edge is
% considered spurious. For example, if you want to consider all traces
% appearing on less than 5 slices, set this value to 5.
%
% OUT
% cleaned_outers: 1xn_slice cell array containing the edge tracings from
% the automated tracing pipeline, with those not appearing on enough slices
% removed
%
%
% R. A. Manzuk 12/21/2020
    %% begin the function
    % get some initial info
    total_edges = max(cellfun(@numel,initial_outers));
    not_empty_ind = zeros(size(initial_outers,2),total_edges);
    
    %figure out which edges are not present in which slices
    for i = 1:numel(initial_outers)
        not_empty_ind(i,1:numel(initial_outers{i})) = ~cellfun(@isempty,initial_outers{i});
    end
    
    %number of slices appearing is just the sum of our not_empty matrix
    not_enough_slices = sum(not_empty_ind) <= nslice_cutoff;
    
    % apply the cutoff
    for i = 1:numel(initial_outers)
        cleaned_outers{i} = initial_outers{i}(~not_enough_slices(1:numel(initial_outers{i})));
    end
end