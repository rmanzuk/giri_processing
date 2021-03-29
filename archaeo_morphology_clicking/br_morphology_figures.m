% Script containing the code to produce final figures for brancing
% morphology data
%
% Ryan A. Manzuk03/23/2021
%% Any 3d data 2d projection
% change these for the data you want to use
data_3d = sm_outer_3d_rotated; 
data_pixel_scale = sm_um_pixel/10000;
center_statistics = sm_outer_center_stats;

figure();
for i = 1:numel(data_3d)
    outline = boundary(data_3d{i}(:,1) .* data_pixel_scale, data_3d{i}(:,3) .* data_pixel_scale);
    pgon = polyshape(data_3d{i}(outline,1) .* data_pixel_scale, data_3d{i}(outline,3) .* data_pixel_scale);
    plot(pgon, 'FaceAlpha', 0.08, 'EdgeAlpha', 0.2, 'FaceColor',[0,0,0]);
    hold on
    plot(center_statistics.spline{i}(:,1) .* data_pixel_scale, center_statistics.spline{i}(:,3) .* data_pixel_scale)
end
xlabel('Modern geographic azimuth')
ylabel('Bedding corrected vertical [cm]')

%% The branch angle histograms

% good to have colors
c_ord  = get(gca,'ColorOrder');

% angles over which we want histograms
angles = [0:180];

% set up some pdfs for other systems
ants_open_pdf = pdf('Normal',angles,43,4);
% ants_forest_pdf = pdf('Normal',angles,65,2.35);
isolated_neurons_pdf = pdf('Normal',angles,98,10);
coral_low_flow = pdf('Normal',angles,90.9,21.9);
%coral_middle_pdf = pdf('Normal',angles,86.2,16.9);
%coral_deep_pdf = pdf('Normal',angles,89.4,13.6);
seepage_channel_pdf = pdf('Normal',angles,72,22.4);

% specify the branch angles for archaeos and corals
% which branch length you would like to use[0.1,0.5,1,2,3] in centimeters
% choose the index
length_used = 4;
sm_angles = sm_br_angles(:,:,length_used);
cc297_angles = cc297_br_angles(:,:,length_used);
caroliana_angles = caroliana_brangles(:,:,length_used);
cytherea_angles = cytherea_brangles(:,:,length_used);
loripes_angles = loripes_brangles(:,:,length_used);
millepora_angles = millepora_brangles(:,:,length_used);

% compile all archaeos and corals into single vectors
archaeo_angles = [unique(sm_angles(sm_angles~=0 & sm_angles<90)); unique(cc297_angles(cc297_angles~=0 & cc297_angles<90))];
coral_angles = [unique(caroliana_angles(caroliana_angles~=0 & caroliana_angles<90));unique(cytherea_angles(cytherea_angles~=0 & cytherea_angles<90));...
    unique(loripes_angles(loripes_angles~=0 & loripes_angles<90));unique(millepora_angles(millepora_angles~=0 & millepora_angles<90))];

gray_level = 1;

figure();
subplot(2,1,1)
area(angles,isolated_neurons_pdf,'DisplayName','Isolated neurons (space filling)','FaceColor',[gray_level,gray_level,gray_level])
hold on
area(angles,coral_low_flow,'DisplayName','Isolated scleractinian (space filling)','FaceColor',[gray_level,gray_level,gray_level],'LineStyle',':')
hold on
area(angles,seepage_channel_pdf,'DisplayName','Seepage channels (diffusion)','FaceColor',[gray_level,gray_level,gray_level],'LineStyle','--')
hold on
area(angles,ants_open_pdf,'DisplayName','Ant trails - open (no material consideration)','FaceColor',[gray_level,gray_level,gray_level],'LineStyle','-.')
hold on
histogram(archaeo_angles,15,'Normalization','pdf', 'DisplayName','Archaeocyathids','FaceColor',c_ord(1,:))
hold on
histogram(coral_angles,15,'Normalization','pdf', 'DisplayName','Reef-building corals','FaceColor',c_ord(2,:))
xlim([0,180])
xlabel('Branching Angle [degrees]')
ylabel('Probability')
legend

all_brangle_data = [archaeo_angles;coral_angles];
labels = [ones(size(unique(sm_angles(sm_angles~=0 & sm_angles<90))));...
    2*ones(size(unique(cc297_angles(cc297_angles~=0 & cc297_angles<90))));...
    4*ones(size(unique(caroliana_angles(caroliana_angles~=0 & caroliana_angles<90))));...
    5*ones(size(unique(cytherea_angles(cytherea_angles~=0 & cytherea_angles<90))));...
    6*ones(size(unique(loripes_angles(loripes_angles~=0 & loripes_angles<90))));...
    7*ones(size(unique(millepora_angles(millepora_angles~=0 & millepora_angles<90))))];
    

subplot(2,1,2)
boxplot(all_brangle_data, labels,'orientation','horizontal')
xlabel('Branching Angle [degrees]')
%% surface area and volume and branch thickness and spacing
sm_vr = (mean(sm_nn_dists(:),'omitnan'))/mean(sm_thicks_encountered(:),'omitnan');
cc297_vr = (mean(cc297_nn_dists(:),'omitnan'))/mean(cc297_thicks_encountered(:),'omitnan');
caroliana_vr = (mean(sm_nn_dists(:),'omitnan'))/mean(caroliana_thicks_encountered(:),'omitnan');
cytherea_vr = (mean(cytherea_nn_dists(:),'omitnan'))/mean(cytherea_thicks_encountered(:),'omitnan');
loripes_vr = (mean(loripes_nn_dists(:),'omitnan'))/mean(loripes_thicks_encountered(:),'omitnan');
millepora_vr = (mean(millepora_nn_dists(:),'omitnan'))/mean(millepora_thicks_encountered(:),'omitnan');

max_number = 2000;
[xx,yy] = meshgrid([1:max_number],[1:max_number]);
slopes = (flip(yy)./xx)/max_number;
figure();
subplot(2,1,1)
pcolor(xx,flip(yy),log2(slopes))
shading flat
cb = colorbar;
colormap(cmocean('rain'))
hold on
scatter(sm_enclosing_volume,sum(sm_surface_area),80,'filled','d', 'MarkerEdgeColor',[1,1,1])
hold on
scatter(cc297_enclosing_volume,sum(cc297_surface_area),80,'filled','d','MarkerEdgeColor',[1,1,1])
hold on
scatter(caroliana_enclosing_volume,sum(caroliana_surface_area), 80,[1,1,1],'filled')
hold on
scatter(cytherea_enclosing_volume,sum(cytherea_surface_area), 80,[1,1,1],'filled')
hold on
scatter(loripes_enclosing_volume,sum(loripes_surface_area), 80,[1,1,1],'filled')
hold on
scatter(millepora_enclosing_volume,sum(millepora_surface_area), 80,[1,1,1],'filled')
xlim([10, 2000])
ylim([10, 2000])
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('enclosing volume cm^{3}')
ylabel('surface area cm^{2}')


max_number = 1.5;
[xx,yy] = meshgrid([0.6:0.01:max_number],[0.01:0.01:0.5]);
slopes = (flip(yy)./xx)/max_number;
subplot(2,1,2)
pcolor(xx,flip(yy),log2(slopes))
shading flat
cb = colorbar;
colormap(cmocean('rain'))
hold on
scatter(mean(sm_nn_dists(:),'omitnan')*sm_um_pixel/1e4,mean(sm_thicks_encountered(:)*sm_um_pixel/1e4,'omitnan'),80,'filled','d', 'MarkerEdgeColor',[1,1,1])
hold on
scatter(mean(cc297_nn_dists(:),'omitnan')*cc297_um_pixel/1e4,mean(cc297_thicks_encountered(:)*cc297_um_pixel/1e4,'omitnan'),80,'filled','d', 'MarkerEdgeColor',[1,1,1])
hold on
scatter(mean(caroliana_nn_dists(:),'omitnan')*caroliana_scale/1e4,mean(caroliana_thicks_encountered(:)*caroliana_scale/1e4,'omitnan'),80,[1,1,1],'filled')
hold on
scatter(mean(cytherea_nn_dists(:),'omitnan')*cytherea_scale/1e4,mean(cytherea_thicks_encountered(:)*cytherea_scale/1e4,'omitnan'),80,[1,1,1],'filled')
hold on
scatter(mean(loripes_nn_dists(:),'omitnan')*loripes_scale/1e4,mean(loripes_thicks_encountered(:)*loripes_scale/1e4,'omitnan'),80,[1,1,1],'filled')
hold on
scatter(mean(millepora_nn_dists(:),'omitnan')*millepora_scale/1e4,mean(millepora_thicks_encountered(:)*millepora_scale/1e4,'omitnan'),80,[1,1,1],'filled')
xlabel('mean branch spacing cm')
ylabel('mean branch radius cm')