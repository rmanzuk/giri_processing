function [surface_area] = std_sa_encvol(outers,scale_ratio, pixel_res, n_cubes, cube_size, branched_flags)
% This function takes outlines of branches that are still separated into
% slices and goes through slice by slice to calcluate surface area and
% volume assuming the branch is cylindrical over that unit. To make for a
% standard unit of measure between samples, the function identifies
% random, cubic region of a specified size and measures SA and Vol in the
% number of regions specified.
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
% pixel_res: unit size of each pixel
%
% n_cubes: Double with the number of randomly placed cubes you would like to sample
%
% cube_size: Double with the size (in original units, pre-unit conversion) of each
% dimension of the cube.
%
% branched_flags: optional, only input if using this for hand clicked
% archaeos. n_branches x n_branches array with the slice numbers for the
% branch points between branches. Output from process_branched_network.m
%
% OUT
% surface_area: row vector of length n_cubes with the total surface 
% area of all branches in each cube. The mean of this vector would be the 
% mean surface area measured.

% R. A. Manzuk 09/20/2021
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
    
    % empty array to get the final surface areas
    surface_area = zeros(1,n_cubes);

    % how thick is each slice in our original units?
    slice_thickness = scale_ratio * pixel_res;

    % how many_slices do we need to make that dimension of the cube?
    slices_needed = ceil(cube_size/slice_thickness);

    % if we don't have enough slices, throw an error
    if slices_needed > n_slices
        error('Error. \nMax cube dimension for thickness is %.3f',(slice_thickness*n_slices))
    end
    
    % and how big is the model in the x-y plane?
    % set up empty array, and just go through the cell array and stack all 
    % points into one 2-column matrix
    all_xy = [];
    for i = 1:numel(outers)
        for j = 1:n_slices
            if ~isempty(outers{i}{j})
                all_xy = [all_xy; outers{i}{j}(:,1:2).*pixel_res];
            end
        end
    end
    
    % then just gather the max, min of those xy coords
    min_xy = min(all_xy);
    max_xy = max(all_xy);
    
    % iterate over the number of cubes
    for q = 1:n_cubes
        % first need to get the random cube
        % now select a random starting point that will let us use that many
        % slices
        start_slice = ceil(rand(1) * (n_slices - slices_needed));
        
        % pick a random xy start that will allow for the cube to fit
        xy_upper_right = (rand(1,2) .* (max_xy-(min_xy+[cube_size,cube_size]))) + min_xy + [cube_size,cube_size];
        min_x = xy_upper_right(1) - cube_size;
        max_x = xy_upper_right(1);
        min_y = xy_upper_right(2) - cube_size;
        max_y = xy_upper_right(2);
        % iterate through all slices
        sa = [];
        for i = start_slice:(start_slice+slices_needed)
            % extract the outline points for each branch in the current slice
            slice_outers = {};
            for j = 1:numel(outers)
                slice_outers(j) = outers{j}(i);
            end
            % if these are archaeo data, the user should input branched flags,
            % which means we need to combine any outlines that have already
            % branched in a smart way that doesn't use their overlapping points in
            % suface area or volume calculations.
            if nargin >5
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
            for j = 1:numel(slice_outers)
                if ~isempty(slice_outers{j})
                    this_outline = slice_outers{j}(:,1:2).*pixel_res;
                    % figure out which portions of branch, if any, are
                    % within the cube
                    in_cube = true(size(this_outline,1),1);
                    in_cube(this_outline(:,1)<min_x) = false;
                    in_cube(this_outline(:,1)>max_x) = false;
                    in_cube(this_outline(:,2)<min_y) = false;
                    in_cube(this_outline(:,2)>max_y) = false;
                    in_cube(isnan(this_outline(:,1))) = false;
                    in_cube(isnan(this_outline(:,2))) = false;
                    % if there's anything left, take the polygon perimeter of the
                    % outline
                    this_outline = this_outline(in_cube,:);
                    if size(this_outline,1)>2
                        warning('off','all')
                        outline_polygon = polyshape(this_outline(:,1),this_outline(:,2));
                        % then translate that perimeter to a surface area and volume of
                        % the segment, assuming a ~cylinder
                        if numel(outline_polygon.Vertices(:,1)) >= 2 && numel(outline_polygon.Vertices(:,2)) >=2
                            perim = arclength(outline_polygon.Vertices(:,1),outline_polygon.Vertices(:,2));
                            if ~isnan(perim)
                                slice_sa = slice_sa + (perim * slice_thickness);
                            end
                        end
                        warning('on','all')
                    end
                end
            end 
            % and update the surface area and volume vectors
            sa = [sa,slice_sa];
        end
        surface_area(q) = sum(sa);
    end
end