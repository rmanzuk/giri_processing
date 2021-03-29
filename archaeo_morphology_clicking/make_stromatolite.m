function [stromat_outlines] = make_stromatolite(radii,height,vert_sampling_rate,circle_sampling_rate,bump_size)
% this function makes the clicking equivalent of a modeled stromatolite.
% Stromatolites are modeled as a vertical stacking of concentric circles
% with varying radii, with some bumps added via sine waves.
%
% IN: 
% radii: row vector with the series of radii for the stromatolite going
% from top to bottom. This vector will be resampled to match the height of
% the stromatolie, so it just needs to capture the relative changes in
% radius. FOr example, a cone could be modeled with radii = [10,8,6,4,2]
%
% height: the desired height of the stromatolite
%
% vert_sampling_rate: number of slices given per unit of height
%
% circle_sampling_rate: number of points given in each slice per unit of
% circumference.
%
% bump_size: size of the bumps
%
% OUT:
%
% stromat_outlines: 1xn_slices cell array where each cell contains the
% x, y coordinates of the outline points for the stromatolite in that
% slice. Multiple stromatolites can then be coalesced into a single cell
% array to compare to branching networks
%
% R. A. Manzuk 03/29/2021
    %% begin the function
    % use sampling rate and height to figure out how many slices we'll have
    n_slices = height * vert_sampling_rate;

    % and resample the radii so be the same size as the n_slices
    rads = imresize(radii,[1,n_slices]);

    % set up an empty cell array to catch the outlines for each slice
    stromat_outlines = cell(1,n_slices);

    % we also want our stromatolite to be a bit bumpy, so we'll apply a sine
    % wave to our circles later. We need to know the amplitude of that sine
    % wave given how big the bumps should be and the vertical position.
    % we'll use a the abs of sine wave that runs the along the vertical of the
    % stromatolite to define the amplitude of the later sine waves
    % samples per unit = vert_sampling_rate
    samples = linspace(0,height,n_slices);
    sine_amp = abs(sin(2*pi*bump_size*samples));


    % allocate the circles given the radius at each height
    for i = 1:n_slices
       % first need to know the circumference to know how many samples to take
       circumf = 2 * pi * rads(i);
       n_points = round(circumf*circle_sampling_rate); 

       % set up theta based upon sampling rate
       theta = linspace(0,2*pi,n_points);

       % apply a sine wave of the proper amplitude to the radius to get bumps
       r = rads(i) + abs(cos(theta*bump_size*sine_amp(i)));

       % now calculate x and y
       x = r .* cos(theta);
       y = r .* sin(theta);

       stromat_outlines{i} = [x',y',zeros(n_points,1)];

    end
end
