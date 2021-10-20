% Script to wrap through functions to select each channel's sharpest image
% and attempt to stack them into a single image
%
%
% R. A. Manzuk 11/20/2020
% last updated 12/08/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set everything up
n_windows = 5;
window_size = 300;
n_zsteps = 100;
z_step_size = 25.4/2;

wavelengths = [470,530,590,625,730,850,940];

p470_dir = '/Users/grinder/Desktop/Data/iq4_experiment_17_/470_output';
p530_dir = '/Users/grinder/Desktop/Data/iq4_experiment_17_/530_output';
p590_dir = '/Users/grinder/Desktop/Data/iq4_experiment_17_/590_output';
p625_dir = '/Users/grinder/Desktop/Data/iq4_experiment_17_/625_output';
p730_dir = '/Users/grinder/Desktop/Data/iq4_experiment_17_/730_output';
p850_dir = '/Users/grinder/Desktop/Data/iq4_experiment_17_/850_output';
p940_dir = '/Users/grinder/Desktop/Data/iq4_experiment_17_/940_output';
%% read all image stacks and get sharpness metrics

[p470_global_sharpness,p470_window_sharpness,edge_coords] = stack_sharpness(p470_dir,5,window_size,n_zsteps,z_step_size);
[p530_global_sharpness,p530_window_sharpness,~] = stack_sharpness(p530_dir,5,window_size,n_zsteps,z_step_size,edge_coords);
[p590_global_sharpness,p590_window_sharpness,~] = stack_sharpness(p590_dir,5,window_size,n_zsteps,z_step_size,edge_coords);
[p625_global_sharpness,p625_window_sharpness,~] = stack_sharpness(p625_dir,5,window_size,n_zsteps,z_step_size,edge_coords);
[p730_global_sharpness,p730_window_sharpness,~] = stack_sharpness(p730_dir,5,window_size,n_zsteps,z_step_size,edge_coords);
[p850_global_sharpness,p850_window_sharpness,~] = stack_sharpness(p850_dir,5,window_size,n_zsteps,z_step_size,edge_coords);
[p940_global_sharpness,p940_window_sharpness,~] = stack_sharpness(p940_dir,5,window_size,n_zsteps,z_step_size,edge_coords);

%% find and read in the sharpest image for each stack...whole vs windows

[~,p470_sharpest_loc] = select_sharpest(p470_dir, p470_global_sharpness, p470_window_sharpness, 'whole',n_zsteps,z_step_size);
[~,p530_sharpest_loc] = select_sharpest(p530_dir, p530_global_sharpness, p530_window_sharpness, 'whole',n_zsteps,z_step_size);
[~,p590_sharpest_loc] = select_sharpest(p470_dir, p590_global_sharpness, p590_window_sharpness, 'whole',n_zsteps,z_step_size);
[~,p625_sharpest_loc] = select_sharpest(p625_dir, p625_global_sharpness, p625_window_sharpness, 'whole',n_zsteps,z_step_size);
[~,p730_sharpest_loc] = select_sharpest(p730_dir, p730_global_sharpness, p730_window_sharpness, 'whole',n_zsteps,z_step_size);
[~,p850_sharpest_loc] = select_sharpest(p850_dir, p850_global_sharpness, p850_window_sharpness, 'whole',n_zsteps,z_step_size);
[~,p940_sharpest_loc] = select_sharpest(p940_dir, p940_global_sharpness, p940_window_sharpness, 'whole',n_zsteps,z_step_size);


%% normalize the sharpness curves

[p470_global_normalized,p470_window_normalized] = normalize_sharpness(p470_global_sharpness,p470_window_sharpness);
[p530_global_normalized,p530_window_normalized] = normalize_sharpness(p530_global_sharpness,p530_window_sharpness);
[p590_global_normalized,p590_window_normalized] = normalize_sharpness(p590_global_sharpness,p590_window_sharpness);
[p625_global_normalized,p625_window_normalized] = normalize_sharpness(p625_global_sharpness,p625_window_sharpness);
[p730_global_normalized,p730_window_normalized] = normalize_sharpness(p730_global_sharpness,p730_window_sharpness);
[p850_global_normalized,p850_window_normalized] = normalize_sharpness(p850_global_sharpness,p850_window_sharpness);
[p940_global_normalized,p940_window_normalized] = normalize_sharpness(p940_global_sharpness,p940_window_sharpness);
%%
subplot(7,1,1)
plot(([1:201]-n_zsteps-1).*z_step_size,p470_global_normalized)
hold on
plot(([1:201]-n_zsteps-1).*z_step_size,p470_window_normalized,'k')
ylabel('Sharpness')
title('p470')
axis tight

subplot(7,1,2)
plot(([1:201]-n_zsteps-1).*z_step_size,p530_global_normalized)
hold on
plot(([1:201]-n_zsteps-1).*z_step_size,p530_window_normalized,'k')
ylabel('Sharpness')
title('p530')
axis tight

subplot(7,1,3)
plot(([1:201]-n_zsteps-1).*z_step_size,p590_global_normalized)
hold on
plot(([1:201]-n_zsteps-1).*z_step_size,p590_window_normalized,'k')
ylabel('Sharpness')
title('p590')
axis tight

subplot(7,1,4)
plot(([1:201]-n_zsteps-1).*z_step_size,p625_global_normalized)
hold on
plot(([1:201]-n_zsteps-1).*z_step_size,p625_window_normalized,'k')
ylabel('Sharpness')
title('p625')
axis tight

subplot(7,1,5)
plot(([1:201]-n_zsteps-1).*z_step_size,p730_global_normalized)
hold on
plot(([1:201]-n_zsteps-1).*z_step_size,p730_window_normalized,'k')
ylabel('Sharpness')
title('p730')
axis tight

subplot(7,1,6)
plot(([1:201]-n_zsteps-1).*z_step_size,p850_global_normalized)
hold on
plot(([1:201]-n_zsteps-1).*z_step_size,p850_window_normalized,'k')
ylabel('Sharpness')
title('p850')
axis tight

subplot(7,1,7)
plot(([1:201]-n_zsteps-1).*z_step_size,p940_global_normalized)
hold on
plot(([1:201]-n_zsteps-1).*z_step_size,p940_window_normalized,'k')
ylabel('Sharpness')
xlabel('z position [um]')
title('p940')
axis tight



%% for the table height vs wavelength plot
sharpness_stretch = 40;


xx = [425:1000];
yy = spline(wavelengths,[p470_sharpest_loc,p530_sharpest_loc,p590_sharpest_loc,p625_sharpest_loc,p730_sharpest_loc,p850_sharpest_loc,p940_sharpest_loc],xx);

plot(wavelengths,[p470_sharpest_loc,p530_sharpest_loc,p590_sharpest_loc,p625_sharpest_loc,p730_sharpest_loc,p850_sharpest_loc,p940_sharpest_loc],'o',xx,yy)
hold on
plot((p470_global_normalized.*sharpness_stretch)+wavelengths(1),([1:201]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((p530_global_normalized.*sharpness_stretch)+wavelengths(2),([1:201]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((p590_global_normalized.*sharpness_stretch)+wavelengths(3),([1:201]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((p625_global_normalized.*sharpness_stretch)+wavelengths(4),([1:201]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((p730_global_normalized.*sharpness_stretch)+wavelengths(5),([1:201]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((p850_global_normalized.*sharpness_stretch)+wavelengths(6),([1:201]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
hold on
plot((p940_global_normalized.*sharpness_stretch)+wavelengths(7),([1:201]-n_zsteps-1).*z_step_size,'k','LineWidth',1)
legend('actual data','spline fit','sharpness curves')
xlabel('wavelength [nm]')
ylabel('focus table height [um]')
title('whole image')
axis tight
ylim([-1000 1000])



%%
green_sharpest_ind = (green_sharpest_loc/z_step_size)+ n_zsteps + 1;

file_pattern = fullfile(green_dir, '*.tif');
tifs = dir(file_pattern);
base_names = natsortfiles({tifs.name});

green_im = imread(fullfile(green_dir, base_names{green_sharpest_ind}));

file_pattern = fullfile(blue_dir, '*.tif');
tifs = dir(file_pattern);
base_names = natsortfiles({tifs.name});

blue_im = imread(fullfile(blue_dir, base_names{green_sharpest_ind}));

%% I made a spreadsheet of all z-stack sharpness: analyze the data.
wavelengths = [470,530,590,625,730,850,940];

[~,sheet_name]=xlsfinfo('/Users/ryan/Dropbox (Princeton)/achromatic_project/z_stack_sharpness.xlsx');
for k=1:numel(sheet_name)
  data{k}=xlsread('/Users/ryan/Dropbox (Princeton)/achromatic_project/z_stack_sharpness.xlsx',sheet_name{k});
end

sharpest_loc = {};
normalized_data = {};
for i = 1:numel(data)
    % extract maxima
    [maxima,sharpest_loc{i}] = max(data{i}(:,4:end));
    % normalize such that all curves go between 0 and 1
    minima = min(data{i}(:,4:end));
    normalized_data{i} = (data{i}(:,4:end) - minima)./(maxima-minima);
end

colors = get(gca,'defaultAxesColorOrder');

figure()
for j = 1:length(wavelengths)
    subplot(ceil(length(wavelengths)/2),2,j)
    for l = 1:numel(data)
        plot((data{l}(:,1)-(data{l}(sharpest_loc{l}(1),1))).*25400,normalized_data{l}(:,(2*j)),'Color',colors(l,:),'DisplayName',sheet_name{l})
        hold on
        plot((data{l}(:,1)-(data{l}(sharpest_loc{l}(1),1))).*25400,normalized_data{l}(:,((2*j)-1)),'Color',colors(l,:),'HandleVisibility','off')
        hold on
    end
    title(num2str(wavelengths(j)))
    xlabel('table position (um)')
    ylabel('sharpness')
    if j == 1
        legend()
    end
end

figure()
for m = 1:numel(data)
    xx = [425:1000];
    global_cols = [1:2:13];
    window_cols = [2:2:14];
    global_sharpest = ((data{m}(sharpest_loc{m}(global_cols)',1)) - data{m}(sharpest_loc{m}(global_cols(2)))).*25400;
    window_sharpest = ((data{m}(sharpest_loc{m}(window_cols)',1)) - data{m}(sharpest_loc{m}(window_cols(2)))).*25400;
    yy_global = spline(wavelengths,global_sharpest,xx);
    yy_window = spline(wavelengths,window_sharpest,xx);
    
    subplot(2,1,1)
    plot(xx,yy_global,'DisplayName',sheet_name{m})
    hold on
    scatter(wavelengths,global_sharpest,'HandleVisibility','off')
    hold on
    title('whole image')
    ylabel('table position (um)')
    xlabel('wavelength (nm)')
    subplot(2,1,2)
    plot(xx,yy_window,'DisplayName',sheet_name{m})
    hold on
    scatter(wavelengths,window_sharpest,'HandleVisibility','off')
    hold on
    title('window mean')
    ylabel('table position (um)')
    xlabel('wavelength (nm)')
    
    
end
legend()

%% those figures are nice to think about changes between samples, but what about 1 figure about chromatic aberration

% we'll use the target because it's flat and easiest to work with (index #
% 5 of the data cell)

% conversion factor between inches and microns
conversion_fac = 25400;

% scale power to accentuate focus differences
scale_power = 200;

% and the limits we'll want for plotting
x_min = 425;
x_max = 1000;
y_min = -500;
y_max = 500; 

% id the columns in the data set that correspond to whole image sharpness
global_cols = [1:2:13];

% then select the table positions (1st column of data) that correspond to
% the sharpest images
[maxima,sharpest_loc] = max(movmean(data{5}(:,4:end),20));
global_sharpest = ((data{5}(sharpest_loc(global_cols)',1)) - data{5}(sharpest_loc(global_cols(2)))).*conversion_fac;

% define the range of x values we want to plot
xx = [x_min:x_max];
% fit a spline for wavelengths and sharpest positions
yy_global = spline(wavelengths,global_sharpest,xx);

% then define each sharpness as a percent of the maximum sharpness
% need the minima
minima = min(movmean(data{5}(:,4:end),20));
pct_max_sharpness = (movmean(data{5}(:,4:end),20)-minima)./(maxima-minima);



% most interested in table heights between y_min and y_max, so define
% that zone

in_window = (data{5}(:,1)- data{5}(sharpest_loc(global_cols(2)))).*conversion_fac > y_min & (data{5}(:,1)- data{5}(sharpest_loc(global_cols(2)))).*conversion_fac < y_max;



% we'll make each wavelengths best sharpness into a patch with alpha
% dictated by sharpness
% the blue patch
patch_x_blue = patch([ones(1,sum(in_window)) .* x_min, ones(1,sum(in_window)) .* x_max],[linspace(y_min,y_max,sum(in_window)), flip(linspace(y_min,y_max,sum(in_window)))],[0,0,1],'HandleVisibility','off');
patch_x_blue.FaceVertexAlphaData = [pct_max_sharpness(in_window,global_cols(1)).^scale_power;flip(pct_max_sharpness(in_window,global_cols(1)).^scale_power)];   
patch_x_blue.FaceAlpha = 'interp' ;
patch_x_blue.EdgeAlpha = 0 ;
hold on

% the green patch
patch_x_green = patch([ones(1,sum(in_window)) .* x_min, ones(1,sum(in_window)) .* x_max],[linspace(y_min,y_max,sum(in_window)), flip(linspace(y_min,y_max,sum(in_window)))],[0,1,0],'HandleVisibility','off');
patch_x_green.FaceVertexAlphaData = [pct_max_sharpness(in_window,global_cols(2)).^scale_power;flip(pct_max_sharpness(in_window,global_cols(2)).^scale_power)];   
patch_x_green.FaceAlpha = 'interp' ;
patch_x_green.EdgeAlpha = 0 ;
hold on

% the yellow patch
patch_x_yellow = patch([ones(1,sum(in_window)) .* x_min, ones(1,sum(in_window)) .* x_max],[linspace(y_min,y_max,sum(in_window)), flip(linspace(y_min,y_max,sum(in_window)))],[1,1,0],'HandleVisibility','off');
patch_x_yellow.FaceVertexAlphaData = [pct_max_sharpness(in_window,global_cols(3)).^scale_power;flip(pct_max_sharpness(in_window,global_cols(3)).^scale_power)];   
patch_x_yellow.FaceAlpha = 'interp' ;  
patch_x_yellow.EdgeAlpha = 0 ;
hold on

% the red patch
patch_x_red = patch([ones(1,sum(in_window)) .* x_min, ones(1,sum(in_window)) .* x_max],[linspace(y_min,y_max,sum(in_window)), flip(linspace(y_min,y_max,sum(in_window)))],[1,0,0],'HandleVisibility','off');
patch_x_red.FaceVertexAlphaData = [pct_max_sharpness(in_window,global_cols(4)).^scale_power;flip(pct_max_sharpness(in_window,global_cols(4)).^scale_power)];   
patch_x_red.FaceAlpha = 'interp' ;  
patch_x_red.EdgeAlpha = 0 ;
hold on

% the red-edge patch
patch_x_red_edge = patch([ones(1,sum(in_window)) .* x_min, ones(1,sum(in_window)) .* x_max],[linspace(y_min,y_max,sum(in_window)), flip(linspace(y_min,y_max,sum(in_window)))],[0.698,0.133,0.133],'HandleVisibility','off');
patch_x_red_edge.FaceVertexAlphaData = [pct_max_sharpness(in_window,global_cols(5)).^scale_power;flip(pct_max_sharpness(in_window,global_cols(5)).^scale_power)];   
patch_x_red_edge.FaceAlpha = 'interp' ;
patch_x_red_edge.EdgeAlpha = 0 ;
hold on

% the nir patch
patch_x_nir = patch([ones(1,sum(in_window)) .* x_min, ones(1,sum(in_window)) .* x_max],[linspace(y_min,y_max,sum(in_window)), flip(linspace(y_min,y_max,sum(in_window)))],[0.5,0.5,0.5],'HandleVisibility','off');
patch_x_nir.FaceVertexAlphaData = [pct_max_sharpness(in_window,global_cols(6)).^scale_power;flip(pct_max_sharpness(in_window,global_cols(6)).^scale_power)];   
patch_x_nir.FaceAlpha = 'interp' ;  
patch_x_nir.EdgeAlpha = 0 ;
hold on

% the ir patch
patch_x_ir = patch([ones(1,sum(in_window)) .* x_min, ones(1,sum(in_window)) .* x_max],[linspace(y_min,y_max,sum(in_window)), flip(linspace(y_min,y_max,sum(in_window)))],[0.1,0.1,0.1]./255,'HandleVisibility','off');
patch_x_ir.FaceVertexAlphaData = [pct_max_sharpness(in_window,global_cols(7)).^scale_power;flip(pct_max_sharpness(in_window,global_cols(7)).^scale_power)];   
patch_x_ir.FaceAlpha = 'interp' ; 
patch_x_ir.EdgeAlpha = 0 ;
hold on



% plot the spline and the scatter points
plot(xx,yy_global, 'k','DisplayName','Spline fit')
hold on
scatter(wavelengths,global_sharpest,40,[0.8,0.8,0.8],'filled','DisplayName','Sharpest image points')
xlabel('wavelength');
ylabel('relative focal distance [\mum]')
axis tight
legend
ylim([-400,300])