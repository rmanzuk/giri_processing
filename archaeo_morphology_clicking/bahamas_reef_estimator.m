% script to get pct reef coverage on bahamas from geyman 2021 figure

% Ryan A. Manzuk 07/09/2021
%% load the image and do some ginput stuff
map_im = imread('geyman_bahamas_figure.png');

imshow(map_im)
% click the left and right limits of the 10km scale bar
[scale_x,~] = ginput(2);

% click the line tracing the barrier reef
[barrier_x,barrier_y] = ginput(2);

% click the perimeter of the sample points
[sample_perim_x,sample_perim_y] = ginput();

% click the perimeter of the whole area
[map_perim_x,map_perim_y] = ginput();
%% make some calculations
% size of scale bar in meters
scale_bar = 10000;
% calculate the per pixel scale
per_pixel_scale = scale_bar/(abs(scale_x(1) - scale_x(2)));
% calculate barrier length, convert to meters, calc an area assuming 10m
% wide
barrier_length = norm([(abs(barrier_x(1) - barrier_x(2))),(abs(barrier_y(1) - barrier_y(2)))]) * per_pixel_scale;
barrier_area = barrier_length * 10;

% make some polygon areas
sample_area = polyarea(sample_perim_x * per_pixel_scale,sample_perim_y * per_pixel_scale);
map_area = polyarea(map_perim_x * per_pixel_scale,map_perim_y * per_pixel_scale); 

% there are 22 blowouts, so what is their percentage area in the sample
% area assuming 100 m^2 per
blowout_area = 22 * 100;
blowout_frac = blowout_area/sample_area;

% same for debris, allowint 10000 m^2 per
debris_area = 22 * 10000;
debris_frac = debris_area/sample_area;

% extrapolate those out to the whole map area
blowout_area_whole = blowout_frac * map_area;
debris_area_whole = debris_frac * map_area;

% percent of map with in situ growth (barrier plus blowouts)
pct_reef = 100 * ((barrier_area + blowout_area_whole)/map_area);

% percent of map with any sign of coral (barrier plus blowouts plus debris)
pct_debris = 100 * ((barrier_area + blowout_area_whole + debris_area_whole)/map_area);

%% process the indexed color tiff of belize
belize_masks = imread('/Users/ryan/Desktop/branch_angle_project/belize_map/indexed_color.tif');
belize_masks = im2double(belize_masks);
belize_masks = rgb2gray(belize_masks);

reef_mask = imbinarize(belize_masks,0.5);
shallow_mask = imbinarize(belize_masks,0.1);

shallow_mask = ~bwareaopen(~shallow_mask,100);

reef_frac = sum(reef_mask)/sum(shallow_mask);