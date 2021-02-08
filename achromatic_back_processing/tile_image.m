function [tiles,centers] = tile_image(image,tile_size)
% this function takes an image and divides it up into tiles of the desired
% size and outputs them as a single 4D array. If your tile size does not
% evenly divide into your image size, the extra bits won't get included in
% the tiles.
%
% IN
%
% image: mxnxn_channels array representing your input image to be tiled
%
% tile_size: single number representing size of square tiles desired or 1x2
% array with rectangular (row,col) size desired
%
% OUT
% 
% tiles: tile_size x tile_size x n_channels x n_tiles array containing the
% output tiles.
%
% tiel_centers: n_tiles x 2 matrix with the row,col indices of the
% locations of each tile's center on the original image. I 
%
% R. A. Manzuk 12/18/2020
    %% begin the function
    if length(tile_size) == 1
        tile_size = [tile_size,tile_size];
    else
    end
    image = im2double(image);   
    tile_index = floor(size(image)./tile_size);

    % which rows and columns will mark the starts and ends of each tile?
    row_starts = [1:tile_size(1):tile_size(1)*tile_index(1)];
    col_starts = [1:tile_size(2):tile_size(2)*tile_index(2)];
    % set up empty arrray to receive tiles
    if ndims(image) == 2
        tiles = zeros(tile_size(1),tile_size(2),tile_index(1)*tile_index(2));
        centers = zeros(tile_index(1)*tile_index(2),2);
        for i = 1:tile_index(1)
            for j = 1:tile_index(2)
                mult_ind = ((i-1)*tile_index(2))+j;
                tiles(:,:,mult_ind) = image(row_starts(i):(row_starts(i)+tile_size(1)-1),col_starts(j):(col_starts(j)+tile_size(2)-1));
                centers(mult_ind,:) = round([row_starts(i)+(tile_size(1)/2),col_starts(j)+(tile_size(2)/2)]);
            end
        end
    else
        tiles = zeros(tile_size(1),tile_size(2), size(image,3),tile_index(1)*tile_index(2));
        centers = zeros(tile_index(1)*tile_index(2),2);
        for i = 1:tile_index(1)
            for j = 1:tile_index(2)
                mult_ind = ((i-1)*tile_index(2))+j;
                tiles(:,:,:,mult_ind) = image(row_starts(i):(row_starts(i)+tile_size(1)-1),col_starts(j):(col_starts(j)+tile_size(2)-1),:);
                centers(mult_ind,:) = round([row_starts(i)+(tile_size(1)/2),col_starts(j)+(tile_size(2)/2)]);
            end
        end
    end
end