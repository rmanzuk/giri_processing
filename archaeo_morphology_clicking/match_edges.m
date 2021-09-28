function [matched_edges] = match_edges(these_edges,previous_edges,distance_thresh)
% This function takes in a tif (uint8 or 16) that is high-contrast and 
% easily segmented as a binary and gets the coordines of the edges of all elements and
% outputs them separately as their own entry in a cell array. Designed to
% mimic hand-clicking function for branching systems but for datasets are high-contrast such that
% outlining can be automated
% IN
%
% these_edges: 1xn_elements cell array where each cell contains the
% coordinates (pixel row and column) for the edges of a single element.
% This is the set of edges that we are looking to match to an existing
% indexing from previous slices.
%
% previous_edges: 1xn_elements cell array where each cell contains the
% coordinates (pixel row and column) for the edges of a single element.
% This is the set of edges already ordered/indexed how we want and will
% look to match with our output
%
% distance_thresh: Euclidean distance at which above which you would say 
% two edges are not the same. 
%
% OUT
% mathed_edges: 1xn_elements cell array where each cell contains the
% coordinates (pixel row and column) for the edges of a single element.
% This set of edges now matches the indexing of the previous edges.
% R. A. Manzuk 12/17/2020
    %% begin the function
    % set up an empty cell array for the matched edges
    matched_edges = cell(1,numel(previous_edges));
    % calculate the centers for this slice and the previous slice
    these_centers = cellfun(@mean,these_edges,'un',0);
    % for the previous one, we need to throw out empty ones
    still_there_ind = ~cellfun(@isempty,previous_edges);
    previous_centers = cellfun(@mean,previous_edges(still_there_ind),'un',0);
    % look at these centers and see if they line up with any of the
    % previous ones
    % easier to work with arrays for now
    these_centers_array = cat(1,these_centers{:});
    previous_centers_array = cat(1,previous_centers{:});
    if numel(previous_centers_array) == 0
        matched_edges = [matched_edges,these_edges];
    else
        for k = 1:size(these_centers_array,1)
            % calculate the distances between this current center and all other
            % previous centers
            distances = zeros(length(still_there_ind),1);
            distances(still_there_ind) = sqrt(sum((previous_centers_array - these_centers_array(k,:)).^2,2));
            % fill in the empty ones with nan
            distances(~still_there_ind) = NaN;
            % find the closest match
            [min_dist,match_ind] = min(distances);
            % if we're close enough, match it
            if min_dist < distance_thresh
                matched_edges{match_ind} = these_edges{k};
            else
                %if we're not close enough, say this is a new element
                current_branches = numel(matched_edges);
                matched_edges{current_branches+1} = these_edges{k};
            end
        end
    end
end