function [resorted_outers] = sort_outline_points(initial_outers)
% this function takes the output of automated edge tracing functions and
% reorders the points so that they go around the edge as opposed to being
% in index order

% IN
% initial_outers: 1xn_slice cell array containing the edge tracings from
% the automated tracing pipeline, could be cleaned up with spurious edges
% removed, and reshaped.
%
% OUT
% resorted_outers: cell array matching dimensions of the input
% containing the edge tracings with points now sorted to go sequentially
% around a circle as opposed to in index order
%
% R. A. Manzuk 12/21/2020
    %% begin the function
    
    % get some initial info
    total_cells = numel(initial_outers);
    total_subcells = max(cellfun(@numel,initial_outers));
    resorted_outers = cell(size(initial_outers));
    
    % loop through and reorder
    for i = 1:total_cells
        for j = 1:total_subcells
                % if there are points there, lets reorder them
            if ~isempty(initial_outers{i}{j})
                points = initial_outers{i}{j};

                % and actually, let's put them in x,y space, not image space
                points_xy = [(points(:,2).*-1),points(:,1)];

                % and then we need to trace the points, starting somewhere and
                % working our way around

                % empty array to accept the resorted points
                new_outers = zeros(size(points_xy));

                % array to tell which points have been sorted as logical
                already_sorted = logical(zeros(1,size(points_xy,1)));

                % we can just keep the first one
                new_outers(1,:) = points_xy(1,:);
                already_sorted(1) = true;

                % and then we need to know the distances between all the points
                distances = pdist2(points_xy,points_xy,'cityblock');

                % and which point what just sorted
                just_sorted = 1;

                for l = 2:length(already_sorted)
                    dist_from_this_point = distances(just_sorted,:);
                    %remove those already sorted (and point itself)
                    dist_from_this_point(already_sorted) = NaN;

                    % which point is closest? update stuff
                    [~,just_sorted] = min(dist_from_this_point);
                    already_sorted(just_sorted) = true;
                    new_outers(l,:) = points_xy(just_sorted,:);
                end 
                resorted_outers{i}{j} = new_outers;
            else
                %if there are no points there, keep it empty
                resorted_outers{i}{j} = [];
            end
        end
    end
end