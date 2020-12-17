function [edge_coords] = get_slice_edges(image,blur_kernel_size,island_size,hole_size)
% This function takes in a tif (uint8 or 16) that is high-contrast and 
% easily segmented as a binary and gets the coordines of the edges of all elements and
% outputs them separately as their own entry in a cell array. Designed to
% mimic hand-clicking function for branching systems but for datasets are high-contrast such that
% outlining can be automated
% IN
%
% image: unit16 or 8 mxnx1 (single-channel) image to be traced
%
% blur_kernel_size: standard deviation of the gaussian to be used for
% smoothing/blurring the image. This will help get rid of small-scale
% features that might complicate things.
%
% island_size: number of pixel threshold for the size of islands to be
% removed in the segmentation (morphological opening)
%
% hole_size: number of pixel threshold for the size of holes to be filled
% in the segmentation (morphological closing)
%
% OUT
% edge_coords: 1xn_elements cell array. each cell then contains the
% coordiantes (pixel row and column) for the edges of a single element.

% R. A. Manzuk 12/17/2020
    %% begin the function
    %make the image blurry, we don't care about the smaller
    blurry_im = imgaussfilt(image,blur_kernel_size);
    % segment just based on threshold
    binary_im = imbinarize(blurry_im);
    % remove islands
    binary_open = bwareaopen(binary_im, island_size,4);
    % close holes
    binary_closed = ~bwareaopen(~binary_open,hole_size,4);
    % detect edges
    edges = edge(binary_closed,'Canny');
    % grab the edge ids and coordinates as connected components
    con_comp = bwconncomp(edges,8);
    % convert the linear indices to sub-indices
    edge_coords = {};
    for j = 1:numel(con_comp.PixelIdxList)
        [edge_coords{j}(:,1),edge_coords{j}(:,2)] = ind2sub(size(image),con_comp.PixelIdxList{j});
    end
end