% script to run through different blocks/functions to analyze different archaeos

% R. A. Manzuk 09/18/2020
%% Load necessary stuff for Stuart's Mill Sample
load('sm_all_inners.mat')
load('sm_all_outers.mat')
block_top_sd = [35,10];
strike_im_heading = 307;
input_folder = '/Users/rmanzuk/Desktop/branch_angle_project/archaeo_clicking_data/grinder_stacks/sm_117_71_downsampled_25';
bedding_sd = [238, 34];
scale_ratio = 6.874/2;
um_pixel = 145.4;

inners = sm_inners;
outers = sm_outers;

%% Load necessary stuff for RLG136a
load('rlg136a_all_inners.mat')
load('rlg136a_all_outers.mat')
scale_ratio = 6.874;
block_top_sd = [194,63];
strike_im_heading = 333;
bedding_sd = [298, 23];
input_folder = '/Users/rmanzuk/Desktop/branch_angle_project/archaeo_clicking_data/grinder_stacks/rlg_136a';
um_pixel = 72.7;

inners = rlg136a_inners;
outers = rlg136a_outers;

%% Load necessary stuff for cc297
load('cc297_all_inners.mat')
load('cc297_all_outers.mat')
scale_ratio = 10;
block_top_sd = [0,90];
strike_im_heading = 270;
bedding_sd = [193, 10];
input_folder = '/Users/rmanzuk/Desktop/branch_angle_project/archaeo_clicking_data/grinder_stacks/cc297/transverse_resample_every10';
um_pixel = 40.4;

inners = cc297_inners;
outers = cc297_outers;

%% Before we do anything, let's densify slices and get the branching points
[inner_dense_slices,outer_dense_slices] = densify_slices(inners,outers,3);
[branched_flags,branching_angles,branch_points_3d] = process_branched_network(inners,outers,scale_ratio,10);

%% Now we can densify in 3d, and take the center lines 
sampling_resolution = 10;
points_here_thresh = 20;
iterate_stop = 4;
sample_freq = 5;
thickness_sampling = 20;
n_iter = 5;

[inner_3d_dense,outer_3d_dense] = densify_3d(inner_dense_slices,outer_dense_slices,scale_ratio,3,0);
[center_points] = iterate_center_lines(outer_3d_dense,sampling_resolution,points_here_thresh,iterate_stop,n_iter,sample_freq,thickness_sampling,15);
%[center_points(38)] = iterate_center_lines(outer_3d_dense(38),3,points_here_thresh,iterate_stop,1,5,10,10);
%[center_points(48)] = iterate_center_lines(outer_3d_dense(48),3,points_here_thresh,iterate_stop,1,5,10,10);
%% and rotate everybody
[centers_rotated,~] = rotate_clicked_data3d(center_points, center_points, block_top_sd, strike_im_heading, bedding_sd,0);
[inner_3d_rotated, outer_3d_rotated] = rotate_clicked_data3d(inner_3d_dense, outer_3d_dense, block_top_sd, strike_im_heading, bedding_sd,1);
[branch_points_rotated] = rotate_branch_points(branch_points_3d, block_top_sd, strike_im_heading, bedding_sd);

%% and get some data
sampling_freq = 5;
thickness_samp = 4;
diff_thresh = 5;
[inner_center_stats,outer_center_stats] = center_line_analysis(inner_3d_rotated,outer_3d_rotated,centers_rotated,sampling_freq,thickness_samp);
[deriv_means,deriv_variances,thicks_encountered] = centers_plane_pass(outer_center_stats,outer_3d_rotated, diff_thresh);
[mean_declinations,mean_slope_runs] = cart2pol(deriv_means(:,1),deriv_means(:,2));
mean_inclinations = atand(deriv_means(:,3)./mean_slope_runs);
%% make some figures
figure(1)
% just the 3d clicking data, rotated
for i = 1:numel(outer_3d_rotated)
plot(outer_3d_rotated{i}(:,1).*um_pixel/10000,outer_3d_rotated{i}(:,3).*um_pixel/10000)
hold on
plot(outer_center_stats.spline{i}(:,1).*um_pixel/10000,outer_center_stats.spline{i}(:,3).*um_pixel/10000);
end
xlabel('Modern geographic azimuth')
ylabel('Bedding-corrected vertical [cm]')

figure(2)
% does thickness depend on inclination
colormap(brewermap(101,'GnBu'))
inc_thick = [];
for i = 1:numel(outer_center_stats.spline)
    inc_thick = [inc_thick;outer_center_stats.inclinations{i},(outer_center_stats.mean_thickness{i}.*um_pixel/1000)'];
end
hist3(inc_thick,'Nbins',[20,20],'CdataMode','auto')
colorbar
view(2)

xlabel('Inclination [degrees]')
ylabel('Mean outer thickness [mm]') 

figure(3)
%2d projections, colorded colored by inclination
subplot(2,2,1)
max_inc = max(cellfun(@max, outer_center_stats.inclinations)); 
min_inc = min(cellfun(@min, outer_center_stats.inclinations)); 
cmap = round(brewermap(101,'GnBu').*255);
for i = 1:numel(outer_center_stats.spline)
    x = outer_center_stats.spline{i}(:,1);
    y = outer_center_stats.spline{i}(:,3);
    incs = outer_center_stats.inclinations{i};
    pct_incs = round(((incs - min_inc)./(max_inc-min_inc)).*100)+1;
    cd = [cmap(pct_incs,:),ones(length(incs),1)];
    p = plot(x,y);
    drawnow
    set(p.Edge,'ColorBinding','interpolated', 'ColorData',uint8(cd'))
    hold on
end
colormap(brewermap(101,'GnBu'))
colorbar
caxis([min_inc, max_inc]);
title('Inclination')

subplot(2,2,2)
%2d projection colored by mean thicknes of the inner skipping 23 
max_thick = max(cellfun(@max, inner_center_stats.mean_thickness)); 
min_thick = min(cellfun(@min, inner_center_stats.mean_thickness)); 
cmap = round(brewermap(101,'GnBu').*255);
for i = [1:numel(outer_center_stats.spline)]
    x = inner_center_stats.spline{i}(:,1);
    y = inner_center_stats.spline{i}(:,3);
    thicks = inner_center_stats.mean_thickness{i};
    thicks(isnan(thicks)) = 2;
    pct_thicks = round(((thicks - min_thick)./(max_thick-min_thick)).*100)+1;
    cd = [cmap(pct_thicks,:),ones(length(thicks),1)];
    p = plot(x,y);
    drawnow
    set(p.Edge,'ColorBinding','interpolated', 'ColorData',uint8(cd'))
    hold on
end
colormap(brewermap(101,'GnBu'))
colorbar
caxis([min_thick*um_pixel/1000, max_thick*um_pixel/1000]);
title('Inner thickness')

subplot(2,2,3)
%2d projection colored by mean thickness
max_thick = max(cellfun(@max, outer_center_stats.mean_thickness)); 
min_thick = min(cellfun(@min, outer_center_stats.mean_thickness)); 
cmap = round(brewermap(101,'GnBu').*255);
for i = 1:numel(outer_center_stats.spline)
    x = outer_center_stats.spline{i}(:,1);
    y = outer_center_stats.spline{i}(:,3);
    thicks = outer_center_stats.mean_thickness{i};
    pct_thicks = round(((thicks - min_thick)./(max_thick-min_thick)).*100)+1;
    cd = [cmap(pct_thicks,:),ones(length(thicks),1)];
    p = plot(x,y);
    drawnow
    set(p.Edge,'ColorBinding','interpolated', 'ColorData',uint8(cd'))
    hold on
end
colormap(brewermap(101,'GnBu'))
colorbar
caxis([min_thick*um_pixel/1000, max_thick*um_pixel/1000]);
title('Outer thickness')

subplot(2,2,4)
%2d projection colored by mean tube thickness
max_thick = max(cellfun(@max, outer_center_stats.mean_thickness)) - min(cellfun(@min, inner_center_stats.mean_thickness)); 
min_thick =  min(cellfun(@min, outer_center_stats.mean_thickness)) - max(cellfun(@max, inner_center_stats.mean_thickness)); 
cmap = round(brewermap(101,'GnBu').*255);
for i = [1:numel(outer_center_stats.spline)]
    x = outer_center_stats.spline{i}(:,1);
    y = outer_center_stats.spline{i}(:,3);
    thicks = outer_center_stats.mean_thickness{i} - inner_center_stats.mean_thickness{i};
    thicks(isnan(thicks)) = 8;
    pct_thicks = round(((thicks - min_thick)./(max_thick-min_thick)).*100)+1;
    cd = [cmap(pct_thicks,:),ones(length(thicks),1)];
    p = plot(x,y);
    drawnow
    set(p.Edge,'ColorBinding','interpolated', 'ColorData',uint8(cd'))
    hold on
end
colormap(brewermap(101,'GnBu'))
colorbar
caxis([min_thick*um_pixel/1000, max_thick*um_pixel/1000]);
title('Tube thickness')
 
% would be nice to just store which parts of network are past branching
% points
above_branch = {};
for i = [1:numel(outer_center_stats.spline)]
    top_branch = max(branch_points_rotated(i,:,3));
    if top_branch == 0
        above_branch{i} = [1:length(outer_center_stats.spline{i}(:,3))];
    else
        those_above = outer_center_stats.spline{i}(:,3) >= top_branch;
        first_ind = find(those_above,1,'first');
        above_branch{i} = [first_ind:length(those_above)];
    end
end

figure(4)
% how thickness changes as we move away from branch points
for i = 1:numel(outer_center_stats.spline)
    if length(above_branch{i}) > 2
        x_vals = outer_center_stats.spline{i}(above_branch{i},1);
        y_vals = outer_center_stats.spline{i}(above_branch{i},2);
        z_vals = outer_center_stats.spline{i}(above_branch{i},2);
        [~,seg_length] = arclength(x_vals,y_vals,z_vals);
        plot(cumsum(seg_length).*um_pixel/1000,outer_center_stats.mean_thickness{i}(above_branch{i}(2:end)).*um_pixel/1000)
        hold on
    else
    end
end
xlabel('Arclength from branching [mm]')
ylabel('Mean outer thickness [mm]') 
ylim([2,4.5])

figure(5)
% tortuosity
for i = 1:numel(outer_center_stats.spline)
    if length(above_branch{i}) > 2
        [arc,~] = arclength(outer_center_stats.spline{i}(above_branch{i},1),outer_center_stats.spline{i}(above_branch{i},2),outer_center_stats.spline{i}(above_branch{i},3));
        p1 = outer_center_stats.spline{i}(above_branch{i}(1),:);
        p2 = outer_center_stats.spline{i}(above_branch{i}(end),:);
        distance = sqrt((p2(1)-p1(1))^2 + (p2(2)-p1(2))^2 + (p2(3)-p1(3))^2); 
        scatter(distance,arc,'filled');
        hold on
    else
    end
    hold on
    plot([0:200],[0:200])
end
ylabel('Total arclength [mm]')
xlabel('straightline distance [mm]') 
%%
% branchgin histos
c_ord  = get(gca,'ColorOrder');
angles = [0:180];
ants_open_pdf = pdf('Normal',angles,43,4);
% ants_forest_pdf = pdf('Normal',angles,65,2.35);
isolated_neurons_pdf = pdf('Normal',angles,98,10);
coral_shallow_pdf = pdf('Normal',angles,90.9,21.9);
%coral_middle_pdf = pdf('Normal',angles,86.2,16.9);
%coral_deep_pdf = pdf('Normal',angles,89.4,13.6);
seepage_channel_pdf = pdf('Normal',angles,72,22.4);

figure(6)
area(angles,isolated_neurons_pdf,'DisplayName','Isolated neurons (space filling)','FaceColor',c_ord(1,:))
hold on
area(angles,coral_shallow_pdf,'DisplayName','Scleractinian (space filling)','FaceColor',c_ord(2,:))
hold on
area(angles,seepage_channel_pdf,'DisplayName','Seepage channels (diffusion)','FaceColor',c_ord(3,:))
% area(angles,ants_forest_pdf,'DisplayName','Ants - forest (material consideration)','FaceColor',c_ord(3,:))
% hold on
area(angles,ants_open_pdf,'DisplayName','Ant trails - open (no material consideration)','FaceColor',c_ord(4,:))
hold on
histogram(unique(branching_angles(branching_angles~=0 & branching_angles<90)),10,'Normalization','pdf','DisplayName','This sample','FaceColor',c_ord(6,:))
% hold on
% plot(angles,coral_middle_pdf,'DisplayName','Middle scleractinian')
% hold on
% plot(angles,coral_deep_pdf,'DisplayName','Deep scleractinian')
xlim([0,180])
xlabel('Branching Angle [degrees]')
ylabel('Probability')
legend
