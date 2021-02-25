function [surface_area,volume] = sa_and_vol(outers,scale_ratio,branched_flags)
% This function takes outlines of branches that are still separated into
% slices and goes through slice by slice to calcluate surface area and
% volume assuming the branch is cylindrical over that unit.
%
% IN
% outers: 1xn_branches cell array containing the cell arrays for outer circles 
% created during data collection. could be from clicking or automated
% tracing.
%
% scale_ratio: ratio of vertical image separation to pixel-width. For
% example if images are separated by 100 microns and pixels are 20 microns,
% this input is 5.
%
% branched_flags: optional, only input if using this for hand clicked
% archaeos. n_branches x n_branches array with the slice numbers for the
% branch points between branches. Output from process_branched_network.m
%
% OUT
% surface_area: column vector of length n_slices with the total surface 
% area of all branches in each slice. Th sum of this vector would be the 
% total surface area
%
% volume: column vector of length n_slices with the total volume of all
% branches in each slice. Th sum of this vector would be the total volume.

% R. A. Manzuk 02/22/2021
    %%
    % find the number of slices in total
    n_slices = max(cellfun(@numel,outers));

    % make sure all individual branch cell arrays are same size...just in
    % case
    for i = 1:numel(outers)
        if numel(outers{i}) < n_slices
            add_empty = cell(1,n_slices - numel(outers{i}));
            outers{i}(numel(outers{i})+1:n_slices) = add_empty;
        else
            % do nothing
        end
    end

    % iterate through all slices
    surface_area = [];
    volume = [];
    for i = 1:n_slices
        % extract the outline points for each branch in the current slice
        slice_outers = {};
        for j = 1:numel(outers)
            slice_outers(j) = outers{j}(i);
        end

        % if these are archaeo data, the user should input branched flags,
        % which means we need to combine any outlines that have already
        % branched in a smart way that doesn't use their overlapping points in
        % suface area or volume calculations.
        if nargin >2
           % figure out which ones we need to worry about branching for
           have_branched = find(branched_flags <i & branched_flags > 0);
           [row_inds,col_inds] = ind2sub(size(branched_flags),have_branched);

           % and some will have branched but are not in this slice
           not_here = find(cellfun(@isempty,slice_outers));
           % so figure out which ones to remove
           row_inds_remove = intersect(row_inds,not_here);
           col_inds_remove = intersect(col_inds,not_here);
           final_removal = or(ismember(row_inds,row_inds_remove),ismember(col_inds,col_inds_remove));
           % remove the emptes and cobine row and column indices
           combined_inds = [row_inds(~final_removal),col_inds(~final_removal)];

           % now,  actually go through and combine those that need it
           % keep track of ones we've gotten rid of 
           trashed_outers = [];
           for j = 1:size(combined_inds,1)
               to_combine = combined_inds(j,:);
               % check that we've not already trashed either one
               if sum(ismember(to_combine,trashed_outers)) == 0
                   % find the union of the 2 polygons
                   warning('off','all')
                   poly1 = polyshape(slice_outers{to_combine(1)}(:,1),slice_outers{to_combine(1)}(:,2));
                   poly2 = polyshape(slice_outers{to_combine(2)}(:,1),slice_outers{to_combine(2)}(:,2));
                   union_poly = union(poly1,poly2);
                   % and remove holes if there
                   if union_poly.NumHoles > 0 
                       union_poly = rmholes(union_poly);
                   end
                   % and then throw these points in one outline and trash the
                   % other
                   slice_outers{to_combine(1)} = union_poly.Vertices;
                   slice_outers{to_combine(2)} = [];
                   trashed_outers = [trashed_outers,to_combine(2)];
                   warning('on','all')
               end
           end
        end

        % now iterate through all of the outlines, get the area and
        % volume present in the slice
        slice_sa = 0;
        slice_vol = 0;
        for j = 1:numel(slice_outers)
            if ~isempty(slice_outers{j})
                % and take the arclength of the outline
                warning('off','all')
                outline_polygon = polyshape(slice_outers{j}(:,1),slice_outers{j}(:,2));
                % then translate that perimeter to a surface area and vulume of
                % the segment, assuming a ~cylinder
                slice_sa = slice_sa + (perimeter(outline_polygon) * scale_ratio);
                slice_vol = slice_vol + (area(outline_polygon) * scale_ratio);
                warning('on','all')
            end
        end 
        % and update teh surface area and volume vectors
        surface_area(i) = slice_sa;
        volume(i) = slice_vol; 

    end
end