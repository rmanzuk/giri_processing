function [cleaned_outers] = remove_spurious_edges(initial_outers)
% This function goes through and removes edges that are outside the ranges
% specified through user inputs on the x-y plane. It helps to get rid of
% elements from the scan that are parts of the coral housing.
%
% IN
% initial_outers: 1xn_slice cell array containing the edge tracings from
% the automated tracing pipeline
%
%
% OUT
% cleaned_outers: 1xn_slice cell array containing the edge tracings from
% the automated tracing pipeline, with bad edges removed
%
%
% R. A. Manzuk 12/21/2020
% edited
    %% begin the function
    % plot a subset of the edges on just an x-y projection to show user
    % which edges they may not want to include
    figure();
    for i = [1:50:numel(initial_outers)]
        for j = 1:numel(initial_outers{i})
            if ~isempty(initial_outers{i}{j})
                scatter(initial_outers{i}{j}(:,1),initial_outers{i}{j}(:,2))
                hold on
            end
        end
    end
    
    % ask the user to identify the cutoff points
    disp('click the spot in the upper left where only spurious branches are up and to the left');
    [x_thresh_ul, y_thresh_ul] = ginput(1);
    disp('click the spot in the upper right where only spurious branches are up and to the right');
    [x_thresh_ur, y_thresh_ur] = ginput(1);
    disp('click the spot in the lower left where only spurious branches are down and to the left');
    [x_thresh_ll, y_thresh_ll] = ginput(1);
    disp('click the spot in the lower right where only spurious branches are down and to the right');
    [x_thresh_lr, y_thresh_lr] = ginput(1);
   

    % now go through the array and figure out which should be excluded
    % based upon thresholds
    % set up empty id system for bad edges
    to_purge = zeros(numel(initial_outers),max(cellfun(@numel,initial_outers)));
    % and an array to keep track of what is empty so we can get rid of
    % edges that just aren't there at all for some reason
    empty_inds = ones(numel(initial_outers),max(cellfun(@numel,initial_outers)));
    
    for i = 1:numel(initial_outers)
        for j = 1:numel(initial_outers{i})
            if ~isempty(initial_outers{i}{j})
                if any(initial_outers{i}{j}(:,1) < x_thresh_ul & initial_outers{i}{j}(:,2) > y_thresh_ul)
                    to_purge(i,j) = 1;
                elseif any(initial_outers{i}{j}(:,1) > x_thresh_ur & initial_outers{i}{j}(:,2) > y_thresh_ur)
                    to_purge(i,j) = 1;
                elseif any(initial_outers{i}{j}(:,1) < x_thresh_ll & initial_outers{i}{j}(:,2) < y_thresh_ll)
                    to_purge(i,j) = 1;
                elseif any(initial_outers{i}{j}(:,1) > x_thresh_lr & initial_outers{i}{j}(:,2) < y_thresh_lr)
                    to_purge(i,j) = 1;
                end
            else
                empty_inds(i,j) = 0;
            end
        end 
    end
    
    % Any edges identified as needing to be purged will have a column sum
    % greater than 0
    bad_edge = sum(to_purge) > 0;
    % any edges that are completely empty will have column sum of 0 in the
    % empty inds matrix
    empty_edge = sum(empty_inds) == 0;
    
    % combine the two above so we know who to get rid of 
    get_rid_inds = or(bad_edge,empty_edge);
    
    % initiate cleaned outers
    % apply the cutoff
    for i = 1:numel(initial_outers)
        cleaned_outers{i} = initial_outers{i}(~get_rid_inds(1:numel(initial_outers{i})));
    end
end