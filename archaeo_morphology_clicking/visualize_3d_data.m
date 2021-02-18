function [] = visualize_3d_data(input_3d, file_name)
% This function takes the branch, cell array of 3d data for this project
% and visualizes it with plot3. If desired by the user, a filename can be
% input to save a gif of the plot.
%
% IN
% input_3d: 1xn_branches cell array with the 3d data for the 3d model being
% made. could be outlines or centerlines or anything.
%
% file_name: (optional) string of the desired file name without extension.
%
% OUT
%
% R. A. Manzuk 02/18/2021
%%
% plot the stuff in 3d
    figure();
    for i = 1:numel(input_3d)
        if ~isempty(input_3d{i})
            plot3(input_3d{i}(:,1),input_3d{i}(:,2),input_3d{i}(:,3))
            hold on
        end
    end
    % if we've specified an output name, we'll make and save a gif.
    if nargin == 2
        % set initail view
        view_azimuth = 0;
        view_elevation = 90;
        view([view_azimuth,view_elevation]);
        % degree step between frames
        degree_step = 5;

        % set up frame grabbing
        frame_count = 71;
        frame = getframe(gcf);
        [im,map] = rgb2ind(frame.cdata,256,'nodither');
        im(1,1,1,frame_count) = 0;
        k = 1;
        % spin 45°
        for i = 0:-degree_step:-45
            view_azimuth = i;
            ([view_azimuth,view_elevation]);
            frame = getframe(gcf);
            im(:,:,1,k) = rgb2ind(frame.cdata,map,'nodither');
            k = k + 1;
        end
        % tilt down
        for i = 90:-degree_step:15
            view_elevation = i;
            view([view_azimuth,view_elevation]);
            frame = getframe(gcf);
            im(:,:,1,k) = rgb2ind(frame.cdata,map,'nodither');
            k = k + 1;
        end
        % spin left
        for i = view_azimuth:-degree_step:-270
            view_azimuth = i;
            view([view_azimuth,view_elevation]);
            frame = getframe(gcf);
            im(:,:,1,k) = rgb2ind(frame.cdata,map,'nodither');
            k = k + 1;
        end
        % spin right
        for i = view_azimuth:degree_step:0
            view_azimuth = i;
            view([view_azimuth,view_elevation]);
            frame = getframe(gcf);
            im(:,:,1,k) = rgb2ind(frame.cdata,map,'nodither');
            k = k + 1;
        end
        % tilt up to original
        for i = view_elevation:degree_step:90
            view_elevation = i;
            view([view_azimuth,view_elevation]);
            frame = getframe(gcf);
            im(:,:,1,k) = rgb2ind(frame.cdata,map,'nodither');
            k = k + 1;
        end
        % and save the gif
        delay_time = 0.1;
        imwrite(im,map,[file_name '.gif'],'DelayTime',delay_time,'LoopCount',inf)
    end
end