% script to make the figure representing scale of different technologies
%
% R. A. Manzuk 07/16/2021
%% list the technologies, their min/max resolution, min/max fov, 
% all listed in microns, res in length of one dimension, fov in squared
% units
nanosims_min_res = 10/512;
nanosims_max_res = 100/256;
nanosims_min_fov = 10^2;
nanosims_max_fov = 100^2;

nanosims_vertices = [nanosims_min_res nanosims_min_fov;...
nanosims_max_res nanosims_min_fov; nanosims_max_res nanosims_max_fov;...
nanosims_min_res nanosims_max_fov];
nanosims_faces = [1 2 3 4];
nanosims_times = [3600;1800;3600;7200];


giri_min_res = 2;
giri_max_res = 3;
giri_min_fov = 4*150000000;
giri_max_fov = 4*(4*150000000);

giri_vertices = [giri_min_res giri_min_fov;...
giri_max_res giri_min_fov; giri_max_res giri_max_fov;...
giri_min_res giri_max_fov];
giri_faces = [1 2 3 4];
giri_times = [10;10;10;10];

all_faces = [1,2,3,4;5,6,7,8];
all_vertices = [nanosims_vertices;giri_vertices];
all_colors = [nanosims_times;giri_times];
%% make the plot
p = patch('Faces',all_faces,'Vertices',all_vertices,'EdgeColor','black','FaceVertexCData',all_colors,'FaceColor','interp');
colormap(brewermap(100,'Blues'));
a = colorbar;
a.Label.String = 'processing time [s]';
xlabel('resolution [\mum]')
ylabel('field of view [\mum^{2}]')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'xdir', 'reverse')
hold on
stipple([nanosims_min_res,nanosims_max_res],[nanosims_min_fov,nanosims_max_fov],1000)
stipple([giri_min_res,giri_max_res],[giri_min_fov,giri_max_fov],100)




%% function to stipple
function [] = stipple(x_dims,y_dims, rel_density)
    a = log10(x_dims(1));
    b = log10(x_dims(2));
    c = log10(y_dims(1));
    d = log10(y_dims(2));
    log_area = (b-a) * (d-c);
    n_points = round(log_area * rel_density);
    x_coords = 10.^(a + (b-a) * rand(1,n_points));
    y_coords = 10.^(c + (d-c) * rand(1,n_points));
    scatter(x_coords,y_coords,0.1,'k')
end


