function [] = visualize_clicking(outers,input_folder)
% This function goes through the clicking stack and outlines and labels the
% archaeos that are a part of the input 'outer' data. The user then 
% clicks through the stack with a slider. It should help see
% what has already been traced and target future branches for tracing
% 
% IN
% outers: 1xn_archaeos cell array containing the cell arrays for outer circles 
% created during data collection. To set this variable up, I just call 
% outers = {outer1, outer2, ..., outern}
% input_folder: pathname string for the folder containing the stack of images
% for a channel
% pause_time: approximate number of seconds you would like to pause on each
% frame to help slow down and get a better look. Somewhere around 0.2 seems
% to work well.
% startNum: what image # you want to start at. Eg if you only wanna
% view tracing data for a subset of the images in your stack
% endNum: what image # in your stack you want to end at. 
%
% OUT
%
% R. A. Manzuk, 07/28/2020
    %% begin the function
    file_pattern = fullfile(input_folder, '*.tif');
    tifs = dir(file_pattern);
    base_names = natsortfiles({tifs.name});
    n_slices = numel(base_names);

    % make sure all individual archaeo cell arrays are same size...just in
        % case
    for i = 1:numel(outers)
        if numel(outers{i}) < n_slices
            add_empty = cell(1,n_slices - numel(outers{i}));
            outers{i}(numel(outers{i})+1:n_slices) = add_empty;
        else
            % do nothing
        end
    end

    n_archaeos = numel(outers);
    slide_step = 1/(n_slices-1);
    
    f = figure('Visible','off');
    c = uicontrol(f,'Style','slider');
    c.Min = 1;
    c.Max = n_slices;
    c.Value = 1;
    c.SliderStep = [slide_step slide_step];
    c.Position = [270 10 60 20];
    c.Callback = @selection;
    f.Visible='on';

    function selection(src,event)
        val = c.Value;
        i = round(val);
        this_im = imread(fullfile(input_folder, base_names{i}));
        slice_outers = {};
        for j = 1:n_archaeos
            slice_outers(j) = outers{j}(i);
        end
        imshow(this_im,'InitialMagnification',400)
        title("image #" + (i));
        hold on
        for j = 1:numel(slice_outers)
            if ~isempty(slice_outers{j})
                warning('off','all');
                poly = polyshape(slice_outers{j}(:,1),slice_outers{j}(:,2));
                warning('on','all');
                [x,y] = centroid(poly);
                txt = string(j);
                plot(slice_outers{j}(:,1),slice_outers{j}(:,2),'r','LineWidth',2)
                text(x,y,txt,'FontSize',14,'Color','w')
                hold on
            else
                % do nothing
            end
        end
    drawnow
    end
end