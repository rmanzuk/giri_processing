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

[~,sheet_name]=xlsfinfo('/Users/rmanzuk/Dropbox (Princeton)/achromatic_project/z_stack_sharpness.xlsx');
for k=1:numel(sheet_name)
  data{k}=xlsread('/Users/rmanzuk/Dropbox (Princeton)/achromatic_project/z_stack_sharpness.xlsx',sheet_name{k});
end

sharpest_loc = {};
normalized_data = {};
for i = 1:numel(data)
    % extract maxima
    [maxima,sharpest_loc{i}] = max(data{i}(:,4:end));
    % normalize
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