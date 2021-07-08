function [edge_coords] = get_slice_edges2(image,blur_kernel_size,island_size,hole_size,thresh)
% This function takes in a double image that is high-contrast and 
% easily segmented as a binary and gets the coordines of the edges of all elements and
% outputs them separately as their own entry in a cell array. Designed to
% mimic hand-clicking function for branching systems but for datasets are high-contrast such that
% outlining can be automated. This is a small adustment to the protocol
% from get_slice_edges to process the CT scans from Kaandorp.

% IN
%
% image: unit16 or 8 mxnx1 (single-channel) image to be traced
%
% blur_kernel_size: standard deviation of the gaussian to be used for
% smoothing/blurring the binary the image. This will help get rid of small-scale
% features that might complicate things.
%
% island_size: number of pixel threshold for the size of islands to be
% removed in the segmentation (morphological opening)
%
% hole_size: number of pixel threshold for the size of holes to be filled
% in the segmentation (morphological closing)
%
% thresh: optional input of the threshold at which you would like to
% binarize the image. If nothing is 
%
% OUT
% edge_coords: 1xn_elements cell array. each cell then contains the
% coordiantes (pixel row and column) for the edges of a single element.

% R. A. Manzuk 04/30/2021
    %% begin the function
    % segment just based on threshold
    if nargin > 4
        binary_im = imbinarize(image,thresh);
    else
        binary_im = imbinarize(image);
    end
    % close holes
    binary_closed = ~bwareaopen(~binary_im, hole_size,4);
    % remove_islands
    binary_opened = bwareaopen(binary_closed,island_size,4);
    
    % and blur the logical to remove craggly edges
    blurred_binary = imgaussfilt(double(binary_opened),blur_kernel_size);
    
    % rebinaraze
    final_binary = imbinarize(blurred_binary,0.7);
    
    % detect edges
    edges = edge(final_binary,'Canny');
    % grab the edge ids and coordinates as connected components
    con_comp = bwconncomp(edges,8);
    % convert the linear indices to sub-indices
    edge_coords = {};
    for j = 1:numel(con_comp.PixelIdxList)
        [edge_coords{j}(:,1),edge_coords{j}(:,2)] = ind2sub(size(image),con_comp.PixelIdxList{j});
    end
    
    % remove edges that are tiny and must be spurious
    size_cells = cellfun(@numel,edge_coords,'UniformOutput',false);
    idx = cellfun(@(x) x <= 10,size_cells,'UniformOutput',false);

    if sum(cell2mat(idx)) > 0
        edge_coords(cell2mat(idx)) = [];
    end
    
    % remove empty cells
    filled_cells = ~cellfun(@isempty,edge_coords);
    edge_coords = edge_coords(filled_cells);
end