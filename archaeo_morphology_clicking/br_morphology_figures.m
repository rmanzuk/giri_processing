% Script containing the code to produce final figures for brancing
% morphology data
%
% Ryan A. Manzuk03/23/2021
%% Any 3d data 2d projection
% change these for the data you want to use
data_3d = cc297_outer_3d_rotated; 
data_pixel_scale = cc297_um_pixel/10000;
center_statistics = cc297_outer_center_stats;

figure();
for i = 1:numel(data_3d)
    outline = boundary(data_3d{i}(:,1) .* data_pixel_scale, data_3d{i}(:,3) .* data_pixel_scale);
    pgon = polyshape(data_3d{i}(outline,1) .* data_pixel_scale, data_3d{i}(outline,3) .* data_pixel_scale);
    plot(pgon, 'FaceAlpha', 0.08, 'EdgeAlpha', 0.2, 'FaceColor',[0,0,0]);
    hold on
    if ~isempty(center_statistics.spline{i})
        plot(center_statistics.spline{i}(:,1) .* data_pixel_scale, center_statistics.spline{i}(:,3) .* data_pixel_scale)
    end
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
ants_forest_pdf = pdf('Normal',angles,65,2.35);
% isolated_neurons_pdf = pdf('Normal',angles,98,10);
%coral_low_flow = pdf('Normal',angles,90.9,21.9);
%coral_middle_pdf = pdf('Normal',angles,86.2,16.9);
%coral_deep_pdf = pdf('Normal',angles,89.4,13.6);
% seepage_channel_pdf = pdf('Normal',angles,72,22.4);

% load the ginput data for streams and neurons
load('/Users/ryan/Desktop/branch_angle_project/branch_angle_gdd/neurons/group_neuron_angles.mat')
load('/Users/ryan/Desktop/branch_angle_project/branch_angle_gdd/neurons/group_neuron_counts.mat')
load('/Users/ryan/Desktop/branch_angle_project/branch_angle_gdd/neurons/iso_neuron_angles.mat')
load('/Users/ryan/Desktop/branch_angle_project/branch_angle_gdd/neurons/iso_neuron_counts.mat')

load('/Users/ryan/Desktop/branch_angle_project/branch_angle_gdd/stream_networks/arid_stream_angles.mat')
load('/Users/ryan/Desktop/branch_angle_project/branch_angle_gdd/stream_networks/arid_stream_freqs.mat')
load('/Users/ryan/Desktop/branch_angle_project/branch_angle_gdd/stream_networks/humid_stream_angles.mat')
load('/Users/ryan/Desktop/branch_angle_project/branch_angle_gdd/stream_networks/humid_stream_freqs.mat')

% compile all archaeos and corals into single vectors
archaeo_angles = [unique(sm_br_angles(sm_br_angles~=0)); unique(cc297_br_angles(cc297_br_angles~=0)); unique(labrador_br_angles(labrador_br_angles~=0 & labrador_br_angles<90))];
coral_angles = [unique(caroliana_brangles(caroliana_brangles~=0));unique(cytherea_brangles(cytherea_brangles~=0));...
    unique(loripes_brangles(loripes_brangles~=0));unique(millepora_brangles(millepora_brangles~=0));...
    unique(madracis6m_brangles(madracis6m_brangles~=0));unique(madracis15m_brangles(madracis15m_brangles~=0));...
    unique(madracis20m_brangles(madracis20m_brangles~=0))];

% make the pdfs
bins = [5:15:170];
archaeo_assignments = histcounts(archaeo_angles,bins);
coral_assignments = histcounts(coral_angles,bins);


archaeo_probs = archaeo_assignments/sum(archaeo_assignments);
coral_probs = coral_assignments/sum(coral_assignments);

bin_midpoints = (bins(2:end) + bins(1:end-1))/2;

archaeo_spline = spline([0,bin_midpoints],[0,archaeo_probs],angles);
coral_spline = spline([0,bin_midpoints],[0,coral_probs],angles);

iso_neuron_spline = spline(iso_neuron_angles,iso_neuron_counts./(sum(iso_neuron_counts)),[30:160]);
group_neuron_spline = spline(group_neuron_angles,group_neuron_counts./(sum(group_neuron_counts)),[30:160]);

% grab the river data from getraer
river_data = readtable('GetraerMaloof2020_data.csv');
river_angles = river_data.angle;
river_aridity = log(river_data.aridity_index);
arid = river_aridity < -1;
humid = river_aridity >0.75;
bins = [5:15:170];
arid_assignments = histcounts(river_angles(arid),bins);
humid_assignments = histcounts(river_angles(humid),bins);
arid_probs = arid_assignments/sum(arid_assignments);
humid_probs = humid_assignments/sum(humid_assignments);
bin_midpoints = (bins(2:end) + bins(1:end-1))/2;
arid_stream_spline = spline([0,bin_midpoints],[0,arid_probs],angles);
humid_stream_spline = spline([0,bin_midpoints],[0,humid_probs],angles);


figure();
subplot(6,1,1)
area([30:160],iso_neuron_spline,'DisplayName','Isolated neurons (space filling)','FaceColor',c_ord(1,:))
xlim([0 180])
ylim([0,max(iso_neuron_spline)])
subplot(6,1,2)
area(angles,arid_stream_spline,'DisplayName','Arid stream network','FaceColor',c_ord(1,:))
hold on
area(angles,humid_stream_spline,'DisplayName','Humid stream network','FaceColor',c_ord(1,:),'LineStyle',':')
xlim([0 180])
ylim([0,max(arid_stream_spline)])
subplot(6,1,3)
area(angles,ants_forest_pdf,'DisplayName','Ant trails - forrest (material consideration)','FaceColor',c_ord(1,:),'LineStyle',':')
hold on
area(angles,ants_open_pdf,'DisplayName','Ant trails - open (no material consideration)','FaceColor',c_ord(1,:))
xlim([0 180])
ylim([0,max(ants_forest_pdf)])
subplot(6,1,4)
area(angles,coral_spline, 'DisplayName','Reef-building corals','FaceColor',c_ord(4,:))
xlim([0,180])
ylim([0,max(coral_spline)])
subplot(6,1,5)
area(angles,archaeo_spline, 'DisplayName','Archaeocyathids','FaceColor',c_ord(7,:))
xlim([0,180])
ylim([0,max(archaeo_spline)])
xlabel('Branching Angle [degrees]')
ylabel('Probability')


all_brangle_data = [archaeo_angles;coral_angles];
labels = [ones(size(unique(sm_br_angles(sm_br_angles~=0))));...
    2*ones(size(unique(cc297_br_angles(cc297_br_angles~=0))));...
    3*ones(size(unique(labrador_br_angles(labrador_br_angles~=0 & labrador_br_angles<90))));...
    5*ones(size(unique(caroliana_brangles(caroliana_brangles~=0))));...
    6*ones(size(unique(cytherea_brangles(cytherea_brangles~=0))));...
    7*ones(size(unique(loripes_brangles(loripes_brangles~=0))));...
    8*ones(size(unique(millepora_brangles(millepora_brangles~=0))));...
    9*ones(size(unique(madracis6m_brangles(madracis6m_brangles~=0))));...
    10*ones(size(unique(madracis15m_brangles(madracis15m_brangles~=0))));...
    11*ones(size(unique(madracis20m_brangles(madracis20m_brangles~=0))))];
    

subplot(6,1,6)
boxplot(all_brangle_data, labels,'orientation','horizontal','symbol','')
xlabel('Branching Angle [degrees]')
xlim([0,180])
%% surface area and volume and branch thickness and spacing
sm_vr = (mean(sm_nn_dists(:),'omitnan'))/mean(sm_thicks_encountered(:),'omitnan');
cc297_vr = (mean(cc297_nn_dists(:),'omitnan'))/mean(cc297_thicks_encountered(:),'omitnan');
labrador_vr = (mean(labrador_nn_dists(:),'omitnan'))/mean(labrador_thicks_encountered(:),'omitnan');
caroliana_vr = (mean(sm_nn_dists(:),'omitnan'))/mean(caroliana_thicks_encountered(:),'omitnan');
cytherea_vr = (mean(cytherea_nn_dists(:),'omitnan'))/mean(cytherea_thicks_encountered(:),'omitnan');
loripes_vr = (mean(loripes_nn_dists(:),'omitnan'))/mean(loripes_thicks_encountered(:),'omitnan');
millepora_vr = (mean(millepora_nn_dists(:),'omitnan'))/mean(millepora_thicks_encountered(:),'omitnan');


colors = brewermap(13,'PuBuGn');

figure()

xs = 1:10000;
plot(xs,0.1*xs,'Color',colors(5,:))
hold on
plot(xs,0.5*xs,'Color',colors(7,:))
hold on
plot(xs,xs,'Color',colors(9,:))
hold on
plot(xs,2*xs,'Color',colors(11,:))
hold on
plot(xs,10*xs,'Color',colors(13,:))
hold on
scatter(caroliana_enclosing_volume,sum(caroliana_surface_area), 80,'filled','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',c_ord(1,:))
hold on
scatter(cytherea_enclosing_volume,sum(cytherea_surface_area), 80,'filled','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',c_ord(2,:))
hold on
scatter(loripes_enclosing_volume,sum(loripes_surface_area), 80,'filled','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',c_ord(3,:))
hold on
scatter(millepora_enclosing_volume,sum(millepora_surface_area), 80,'filled','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',c_ord(4,:))
hold on
scatter(sm_enclosing_volume,sum(sm_surface_area),80,'filled','d', 'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',c_ord(5,:))
hold on
scatter(cc297_enclosing_volume,sum(cc297_surface_area),80,'filled','d','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',c_ord(6,:))
hold on
scatter(labrador_enclosing_volume,sum(labrador_surface_area),80,'filled','d','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',c_ord(7,:))
hold on
scatter(stromat_enclosing_volume,sum(stromat_surface_area),80,'filled','s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',c_ord(8,:))
xlim([10, 10000])
ylim([10, 2000])
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('enclosing volume [cm^{3}]')
ylabel('surface area [cm^{2}]')
%% weighted least squares
xdata = [median(caroliana_nn_dists(:),'omitnan')*caroliana_scale/1e4, median(cytherea_nn_dists(:),'omitnan')*cytherea_scale/1e4, median(loripes_nn_dists(:),'omitnan')*loripes_scale/1e4,...
    median(millepora_nn_dists(:),'omitnan')*millepora_scale/1e4, median(madracis6m_nn_dists(:),'omitnan')*madracis6m_scale/1e4,...
    median(madracis15m_nn_dists(:),'omitnan')*madracis15m_scale/1e4, median(madracis20m_nn_dists(:),'omitnan')*madracis20m_scale/1e4,...
    median(sm_nn_dists(:),'omitnan')*sm_um_pixel/1e4, median(cc297_nn_dists(:),'omitnan')*cc297_um_pixel/1e4, median(labrador_nn_dists(:),'omitnan')*labrador_um_pixel/1e4];

ydata = [median(caroliana_thicks_encountered(:),'omitnan')*caroliana_scale/1e4, median(cytherea_thicks_encountered(:),'omitnan')*cytherea_scale/1e4, median(loripes_thicks_encountered(:),'omitnan')*loripes_scale/1e4,...
    median(millepora_thicks_encountered(:),'omitnan')*millepora_scale/1e4, median(madracis6m_thicks_encountered(:),'omitnan')*madracis6m_scale/1e4,...
    median(madracis15m_thicks_encountered(:),'omitnan')*madracis15m_scale/1e4, median(madracis20m_thicks_encountered(:),'omitnan')*madracis20m_scale/1e4,...
    median(sm_thicks_encountered(:),'omitnan')*sm_um_pixel/1e4, median(cc297_thicks_encountered(:),'omitnan')*cc297_um_pixel/1e4, median(labrador_thicks_encountered(:),'omitnan')*labrador_um_pixel/1e4];

std_weights_x = [std(mean(caroliana_nn_dists,'omitnan'),'omitnan')*caroliana_scale/1e4, std(mean(cytherea_nn_dists,'omitnan'),'omitnan')*cytherea_scale/1e4,...
    std(mean(loripes_nn_dists,'omitnan'),'omitnan')*loripes_scale/1e4, std(mean(millepora_nn_dists,'omitnan'),'omitnan')*millepora_scale/1e4,...
    std(mean(madracis6m_nn_dists,'omitnan'),'omitnan')*madracis6m_scale/1e4, std(mean(madracis15m_nn_dists,'omitnan'),'omitnan')*madracis15m_scale/1e4,...
    std(mean(madracis20m_nn_dists,'omitnan'),'omitnan')*madracis20m_scale/1e4, std(mean(sm_nn_dists,'omitnan'),'omitnan')*sm_um_pixel/1e4,...
    std(mean(cc297_nn_dists,'omitnan'),'omitnan')*cc297_um_pixel/1e4, std(mean(labrador_nn_dists,'omitnan'),'omitnan')*labrador_um_pixel/1e4];

std_weights_y = [std(mean(caroliana_thicks_encountered,'omitnan'),'omitnan')*caroliana_scale/1e4, std(mean(cytherea_thicks_encountered,'omitnan'),'omitnan')*cytherea_scale/1e4,...
    std(mean(loripes_thicks_encountered,'omitnan'),'omitnan')*loripes_scale/1e4, std(mean(millepora_thicks_encountered,'omitnan'),'omitnan')*millepora_scale/1e4,...
    std(mean(madracis6m_thicks_encountered,'omitnan'),'omitnan')*madracis6m_scale/1e4, std(mean(madracis15m_thicks_encountered,'omitnan'),'omitnan')*madracis15m_scale/1e4,...
    std(mean(madracis20m_thicks_encountered,'omitnan'),'omitnan')*madracis20m_scale/1e4, std(mean(sm_thicks_encountered,'omitnan'),'omitnan')*sm_um_pixel/1e4,...
    std(mean(cc297_thicks_encountered,'omitnan'),'omitnan')*cc297_um_pixel/1e4, std(mean(labrador_thicks_encountered,'omitnan'),'omitnan')*labrador_um_pixel/1e4];
covary = cov(xdata',ydata');
pearson = covary./(std(xdata)*std(ydata));
%%

% load kaandorp data
kaandorp_br_data = readtable('/Users/ryan/Desktop/branch_angle_project/coral_ct_data/kaandorp_scans/filatov2013_coral_data.csv');



low_pct = 40;
high_pct = 60;
figure()
scatter(kaandorp_br_data.br_spacing_mean./5,kaandorp_br_data.db_mean./20)
hold on
scatter(median(caroliana_nn_dists(:),'omitnan')*caroliana_scale/1e4,median(caroliana_thicks_encountered(:)*caroliana_scale/1e4,'omitnan'),80,'filled', 'MarkerEdgeColor',[0,0,0])
dist_prct_l = prctile(caroliana_nn_dists,low_pct,'all')*caroliana_scale/1e4;
dist_prct_h = prctile(caroliana_nn_dists,high_pct,'all')*caroliana_scale/1e4;
rad_prct_l = prctile(caroliana_thicks_encountered,low_pct,'all')*caroliana_scale/1e4;
rad_prct_h = prctile(caroliana_thicks_encountered,high_pct,'all')*caroliana_scale/1e4;
plotEllipses([(dist_prct_l + dist_prct_h)/2,(rad_prct_l + rad_prct_h)/2],...
    [(dist_prct_h - dist_prct_l)/2,(rad_prct_h - rad_prct_l)/2]);

hold on

scatter(median(cytherea_nn_dists(:),'omitnan')*cytherea_scale/1e4,median(cytherea_thicks_encountered(:)*cytherea_scale/1e4,'omitnan'),80,'filled', 'MarkerEdgeColor',[0,0,0])
dist_prct_l = prctile(cytherea_nn_dists,low_pct,'all')*cytherea_scale/1e4;
dist_prct_h = prctile(cytherea_nn_dists,high_pct,'all')*cytherea_scale/1e4;
rad_prct_l = prctile(cytherea_thicks_encountered,low_pct,'all')*cytherea_scale/1e4;
rad_prct_h = prctile(cytherea_thicks_encountered,high_pct,'all')*cytherea_scale/1e4;
plotEllipses([(dist_prct_l + dist_prct_h)/2,(rad_prct_l + rad_prct_h)/2],...
    [(dist_prct_h - dist_prct_l)/2,(rad_prct_h - rad_prct_l)/2]);

hold on

scatter(median(loripes_nn_dists(:),'omitnan')*loripes_scale/1e4,median(loripes_thicks_encountered(:)*loripes_scale/1e4,'omitnan'),80,'filled', 'MarkerEdgeColor',[0,0,0])
dist_prct_l = prctile(loripes_nn_dists,low_pct,'all')*loripes_scale/1e4;
dist_prct_h = prctile(loripes_nn_dists,high_pct,'all')*loripes_scale/1e4;
rad_prct_l = prctile(loripes_thicks_encountered,low_pct,'all')*loripes_scale/1e4;
rad_prct_h = prctile(loripes_thicks_encountered,high_pct,'all')*loripes_scale/1e4;
plotEllipses([(dist_prct_l + dist_prct_h)/2,(rad_prct_l + rad_prct_h)/2],...
    [(dist_prct_h - dist_prct_l)/2,(rad_prct_h - rad_prct_l)/2]);

hold on

scatter(median(millepora_nn_dists(:),'omitnan')*millepora_scale/1e4,median(millepora_thicks_encountered(:)*millepora_scale/1e4,'omitnan'),80,'filled', 'MarkerEdgeColor',[0,0,0])
dist_prct_l = prctile(millepora_nn_dists,low_pct,'all')*millepora_scale/1e4;
dist_prct_h = prctile(millepora_nn_dists,high_pct,'all')*millepora_scale/1e4;
rad_prct_l = prctile(millepora_thicks_encountered,low_pct,'all')*millepora_scale/1e4;
rad_prct_h = prctile(millepora_thicks_encountered,high_pct,'all')*millepora_scale/1e4;
plotEllipses([(dist_prct_l + dist_prct_h)/2,(rad_prct_l + rad_prct_h)/2],...
    [(dist_prct_h - dist_prct_l)/2,(rad_prct_h - rad_prct_l)/2]);

hold on

scatter(median(madracis6m_nn_dists(:),'omitnan')*madracis6m_scale/1e4,median(madracis6m_thicks_encountered(:)*madracis6m_scale/1e4,'omitnan'),80,'filled', 'MarkerEdgeColor',[0,0,0])
dist_prct_l = prctile(madracis6m_nn_dists,low_pct,'all')*madracis6m_scale/1e4;
dist_prct_h = prctile(madracis6m_nn_dists,high_pct,'all')*madracis6m_scale/1e4;
rad_prct_l = prctile(madracis6m_thicks_encountered,low_pct,'all')*madracis6m_scale/1e4;
rad_prct_h = prctile(madracis6m_thicks_encountered,high_pct,'all')*madracis6m_scale/1e4;
plotEllipses([(dist_prct_l + dist_prct_h)/2,(rad_prct_l + rad_prct_h)/2],...
    [(dist_prct_h - dist_prct_l)/2,(rad_prct_h - rad_prct_l)/2]);

hold on

scatter(median(madracis15m_nn_dists(:),'omitnan')*madracis15m_scale/1e4,median(madracis15m_thicks_encountered(:)*madracis15m_scale/1e4,'omitnan'),80,'filled', 'MarkerEdgeColor',[0,0,0])
dist_prct_l = prctile(madracis15m_nn_dists,low_pct,'all')*madracis15m_scale/1e4;
dist_prct_h = prctile(madracis15m_nn_dists,high_pct,'all')*madracis15m_scale/1e4;
rad_prct_l = prctile(madracis15m_thicks_encountered,low_pct,'all')*madracis15m_scale/1e4;
rad_prct_h = prctile(madracis15m_thicks_encountered,high_pct,'all')*madracis15m_scale/1e4;
plotEllipses([(dist_prct_l + dist_prct_h)/2,(rad_prct_l + rad_prct_h)/2],...
    [(dist_prct_h - dist_prct_l)/2,(rad_prct_h - rad_prct_l)/2]);

hold on

scatter(median(madracis20m_nn_dists(:),'omitnan')*madracis20m_scale/1e4,median(madracis20m_thicks_encountered(:)*madracis20m_scale/1e4,'omitnan'),80,'filled', 'MarkerEdgeColor',[0,0,0])
dist_prct_l = prctile(madracis20m_nn_dists,low_pct,'all')*madracis20m_scale/1e4;
dist_prct_h = prctile(madracis20m_nn_dists,high_pct,'all')*madracis20m_scale/1e4;
rad_prct_l = prctile(madracis20m_thicks_encountered,low_pct,'all')*madracis20m_scale/1e4;
rad_prct_h = prctile(madracis20m_thicks_encountered,high_pct,'all')*madracis20m_scale/1e4;
plotEllipses([(dist_prct_l + dist_prct_h)/2,(rad_prct_l + rad_prct_h)/2],...
    [(dist_prct_h - dist_prct_l)/2,(rad_prct_h - rad_prct_l)/2]);

hold on

scatter(median(sm_nn_dists(:),'omitnan')*sm_um_pixel/1e4,median(sm_thicks_encountered(:)*sm_um_pixel/1e4,'omitnan'),80,'filled','d', 'MarkerEdgeColor',[0,0,0])
dist_prct_l = prctile(sm_nn_dists,low_pct,'all')*sm_um_pixel/1e4;
dist_prct_h = prctile(sm_nn_dists,high_pct,'all')*sm_um_pixel/1e4;
rad_prct_l = prctile(sm_thicks_encountered,low_pct,'all')*sm_um_pixel/1e4;
rad_prct_h = prctile(sm_thicks_encountered,high_pct,'all')*sm_um_pixel/1e4;
plotEllipses([(dist_prct_l + dist_prct_h)/2,(rad_prct_l + rad_prct_h)/2],...
    [(dist_prct_h - dist_prct_l)/2,(rad_prct_h - rad_prct_l)/2]);

hold on

scatter(median(cc297_nn_dists(:),'omitnan')*cc297_um_pixel/1e4,median(cc297_thicks_encountered(:)*cc297_um_pixel/1e4,'omitnan'),80,'filled','d', 'MarkerEdgeColor',[0,0,0])
dist_prct_l = prctile(cc297_nn_dists,low_pct,'all')*cc297_um_pixel/1e4;
dist_prct_h = prctile(cc297_nn_dists,high_pct,'all')*cc297_um_pixel/1e4;
rad_prct_l = prctile(cc297_thicks_encountered,low_pct,'all')*cc297_um_pixel/1e4;
rad_prct_h = prctile(cc297_thicks_encountered,high_pct,'all')*cc297_um_pixel/1e4;
plotEllipses([(dist_prct_l + dist_prct_h)/2,(rad_prct_l + rad_prct_h)/2],...
    [(dist_prct_h - dist_prct_l)/2,(rad_prct_h - rad_prct_l)/2]);

hold on

scatter(median(labrador_nn_dists(:),'omitnan')*labrador_um_pixel/1e4,median(labrador_thicks_encountered(:)*labrador_um_pixel/1e4,'omitnan'),80,'filled','d', 'MarkerEdgeColor',[0,0,0])
dist_prct_l = prctile(labrador_nn_dists,low_pct,'all')*labrador_um_pixel/1e4;
dist_prct_h = prctile(labrador_nn_dists,high_pct,'all')*labrador_um_pixel/1e4;
rad_prct_l = prctile(labrador_thicks_encountered,low_pct,'all')*labrador_um_pixel/1e4;
rad_prct_h = prctile(labrador_thicks_encountered,high_pct,'all')*labrador_um_pixel/1e4;
plotEllipses([(dist_prct_l + dist_prct_h)/2,(rad_prct_l + rad_prct_h)/2],...
    [(dist_prct_h - dist_prct_l)/2,(rad_prct_h - rad_prct_l)/2]);

xlabel('branch spacing [cm]')
ylabel('branch radius [cm]')
%%
plot(unique_ages,sa_vol_bins)
hold on
scatter(520,sum(sm_surface_area)/sm_enclosing_volume)
hold on
scatter(520,sum(cc297_surface_area)/cc297_enclosing_volume)
hold on
scatter(520,sum(labrador_surface_area)/labrador_enclosing_volume)
hold on
scatter(400,sum(millepora_surface_area)/millepora_enclosing_volume)
hold on
scatter(400,sum(loripes_surface_area)/loripes_enclosing_volume)
hold on
scatter(400,sum(cytherea_surface_area)/cytherea_enclosing_volume)
hold on
scatter(400,sum(caroliana_surface_area)/caroliana_enclosing_volume)
hold on
scatter(400,sum(madracis20m_surface_area)/madracis20m_enclosing_volume)
hold on
scatter(400,sum(madracis15m_surface_area)/madracis15m_enclosing_volume)
hold on
scatter(400,sum(madracis6m_surface_area)/madracis6m_enclosing_volume)
xlabel('Time [Ma]')
ylabel('Surface area / enclosing volume [cm^{-1}]')
set(gca, 'xdir', 'reverse')
%%
figure()
for i = 1:size(c_ord,1)
   scatter(1,i,'filled')
   hold on
end