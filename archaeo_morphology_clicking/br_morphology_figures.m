% Script containing the code to produce final figures for brancing
% morphology data
%
% Ryan A. Manzuk03/23/2021
%% Any 3d data 2d projection
% change these for the data you want to use
data_3d = madracis6m_3d; 
data_pixel_scale = 250/1e4;
center_statistics = madracis6m_center_stats;

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

%% make a bunch of box & whisker plots
figure()
subplot(2,1,1)

ant_open_pseudo = randn(1000,1) * 4 + 43;
ant_forest_pseudo = randn(1000,1) * 2.35 + 65;

all_brangle_data = [archaeo_angles;coral_angles;river_angles(arid);river_angles(humid);...
    ant_open_pseudo;ant_forest_pseudo];
labels = [ones(size(unique(sm_br_angles(sm_br_angles~=0))));...
    2*ones(size(unique(cc297_br_angles(cc297_br_angles~=0))));...
    3*ones(size(unique(labrador_br_angles(labrador_br_angles~=0 & labrador_br_angles<90))));...
    5*ones(size(unique(caroliana_brangles(caroliana_brangles~=0))));...
    6*ones(size(unique(cytherea_brangles(cytherea_brangles~=0))));...
    7*ones(size(unique(loripes_brangles(loripes_brangles~=0))));...
    8*ones(size(unique(millepora_brangles(millepora_brangles~=0))));...
    9*ones(size(unique(madracis6m_brangles(madracis6m_brangles~=0))));...
    10*ones(size(unique(madracis15m_brangles(madracis15m_brangles~=0))));...
    11*ones(size(unique(madracis20m_brangles(madracis20m_brangles~=0))));...
    13*ones(sum(arid),1);14*ones(sum(humid),1);15*ones(numel(ant_open_pseudo),1);16*ones(numel(ant_forest_pseudo),1)];
boxplot(all_brangle_data, labels,'orientation','horizontal','symbol','')
xlabel('Branching Angle [degrees]')
xlim([0,180])

subplot(2,3,4)

sm_head_angles = heading_angles(sm_outer_center_stats);
cc297_head_angles = heading_angles(cc297_outer_center_stats);
labrador_head_angles = heading_angles(labrador_outer_center_stats);
caroliana_head_angles = heading_angles(caroliana_center_stats);
cytherea_head_angles = heading_angles(cytherea_center_stats);
loripes_head_angles = heading_angles(loripes_center_stats);
millepora_head_angles = heading_angles(millepora_center_stats);
madracis6m_head_angles = heading_angles(madracis6m_center_stats);
madracis15m_head_angles = heading_angles(madracis15m_center_stats);
madracis20m_head_angles = heading_angles(madracis20m_center_stats);

% all_heading_data = [sm_head_angles';cc297_head_angles';labrador_head_angles';caroliana_head_angles'...
%     ;cytherea_head_angles';loripes_head_angles';millepora_head_angles';madracis6m_head_angles'...
%     ;madracis15m_head_angles';madracis20m_head_angles'];
% labels = [ones(numel(sm_head_angles),1);...
%     2*ones(numel(cc297_head_angles),1);...
%     3*ones(numel(labrador_head_angles),1);...
%     5*ones(numel(caroliana_head_angles),1);...
%     6*ones(numel(cytherea_head_angles),1);...
%     7*ones(numel(loripes_head_angles),1);...
%     8*ones(numel(millepora_head_angles),1);...
%     9*ones(numel(madracis6m_head_angles),1);...
%     10*ones(numel(madracis15m_head_angles),1);...
%     11*ones(numel(madracis20m_head_angles),1)];
all_heading_data = [(mean(sm_deriv_variances,2,'omitnan'));(mean(cc297_deriv_variances,2,'omitnan'));...
    (mean(labrador_deriv_variances,2,'omitnan'));(mean(caroliana_deriv_variances,2,'omitnan'))...
    ;(mean(cytherea_deriv_variances,2,'omitnan'));(mean(loripes_deriv_variances,2,'omitnan'));...
    (mean(millepora_deriv_variances,2,'omitnan'));(mean(madracis6m_deriv_variances,2,'omitnan'))...
    ;(mean(madracis15m_deriv_variances,2,'omitnan'));(mean(madracis20m_deriv_variances,2,'omitnan'))];
labels = [ones(numel(mean(sm_deriv_variances,2,'omitnan')),1);...
    2*ones(numel(mean(cc297_deriv_variances,2,'omitnan')),1);...
    3*ones(numel(mean(labrador_deriv_variances,2,'omitnan')),1);...
    5*ones(numel(mean(caroliana_deriv_variances,2,'omitnan')),1);...
    6*ones(numel(mean(cytherea_deriv_variances,2,'omitnan')),1);...
    7*ones(numel(mean(loripes_deriv_variances,2,'omitnan')),1);...
    8*ones(numel(mean(millepora_deriv_variances,2,'omitnan')),1);...
    9*ones(numel(mean(madracis6m_deriv_variances,2,'omitnan')),1);...
    10*ones(numel(mean(madracis15m_deriv_variances,2,'omitnan')),1);...
    11*ones(numel(mean(madracis20m_deriv_variances,2,'omitnan')),1)];
boxplot(all_heading_data, labels,'symbol','')
ylabel('Heading variation')
ylim([0,3])

subplot(2,3,5)

all_radii_data = [mean(sm_thicks_encountered,'omitnan')'*sm_um_pixel/1e4;mean(cc297_thicks_encountered,'omitnan')'*cc297_um_pixel/1e4;...
    mean(labrador_thicks_encountered,'omitnan')'*labrador_um_pixel/1e4;mean(caroliana_thicks_encountered,'omitnan')'*caroliana_scale/1e4;...
    mean(cytherea_thicks_encountered,'omitnan')'*cytherea_scale/1e4;mean(loripes_thicks_encountered,'omitnan')'*loripes_scale/1e4;...
    mean(millepora_thicks_encountered,'omitnan')'*millepora_scale/1e4;mean(madracis6m_thicks_encountered,'omitnan')'*madracis6m_scale/1e4;...
    mean(madracis15m_thicks_encountered,'omitnan')'*madracis15m_scale/1e4;mean(madracis20m_thicks_encountered,'omitnan')'*madracis20m_scale/1e4];
labels = [ones(numel(mean(sm_thicks_encountered,'omitnan')),1);...
    2*ones(numel(mean(cc297_thicks_encountered,'omitnan')),1);...
    3*ones(numel(mean(labrador_thicks_encountered,'omitnan')),1);...
    5*ones(numel(mean(caroliana_thicks_encountered,'omitnan')),1);...
    6*ones(numel(mean(cytherea_thicks_encountered,'omitnan')),1);...
    7*ones(numel(mean(loripes_thicks_encountered,'omitnan')),1);...
    8*ones(numel(mean(millepora_thicks_encountered,'omitnan')),1);...
    9*ones(numel(mean(madracis6m_thicks_encountered,'omitnan')),1);...
    10*ones(numel(mean(madracis15m_thicks_encountered,'omitnan')),1);...
    11*ones(numel(mean(madracis20m_thicks_encountered,'omitnan')),1)];
boxplot(all_radii_data, labels,'symbol','')
ylabel('Branch radius [cm]')
ylim([0,1])

subplot(2,3,6)

all_spacing_data = [mean(sm_nn_dists,'omitnan')'*sm_um_pixel/1e4;mean(cc297_nn_dists,'omitnan')'*cc297_um_pixel/1e4;...
    mean(labrador_nn_dists,'omitnan')'*labrador_um_pixel/1e4;mean(caroliana_nn_dists,'omitnan')'*caroliana_scale/1e4;...
    mean(cytherea_nn_dists,'omitnan')'*cytherea_scale/1e4;mean(loripes_nn_dists,'omitnan')'*loripes_scale/1e4;...
    mean(millepora_nn_dists,'omitnan')'*millepora_scale/1e4;mean(madracis6m_nn_dists,'omitnan')'*madracis6m_scale/1e4;...
    mean(madracis15m_nn_dists,'omitnan')'*madracis15m_scale/1e4;mean(madracis20m_nn_dists,'omitnan')'*madracis20m_scale/1e4];
labels = [ones(numel(mean(sm_nn_dists,'omitnan')),1);...
    2*ones(numel(mean(cc297_nn_dists,'omitnan')),1);...
    3*ones(numel(mean(labrador_nn_dists,'omitnan')),1);...
    5*ones(numel(mean(caroliana_nn_dists,'omitnan')),1);...
    6*ones(numel(mean(cytherea_nn_dists,'omitnan')),1);...
    7*ones(numel(mean(loripes_nn_dists,'omitnan')),1);...
    8*ones(numel(mean(millepora_nn_dists,'omitnan')),1);...
    9*ones(numel(mean(madracis6m_nn_dists,'omitnan')),1);...
    10*ones(numel(mean(madracis15m_nn_dists,'omitnan')),1);...
    11*ones(numel(mean(madracis20m_nn_dists,'omitnan')),1)];
boxplot(all_spacing_data, labels,'symbol','')
ylabel('Branch spacing [cm]')
ylim([0,3])
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

lower_pct = 25;
upper_pct = 75;
plot(unique_ages,sa_vol_bins*10000/cube_size^3)
hold on
scatter(520,mean(sm_surf_areas)*10000/cube_size^3,'d')
errorbar(520,mean(sm_surf_areas)*10000/cube_size^3,...
    mean(sm_surf_areas)*10000/cube_size^3-prctile(sm_surf_areas,lower_pct)*10000/cube_size^3,...
    prctile(sm_surf_areas,upper_pct)*10000/cube_size^3-mean(sm_surf_areas)*10000/cube_size^3)
hold on
scatter(520,mean(cc297_surf_areas)*10000/cube_size^3)
errorbar(520,mean(cc297_surf_areas)*10000/cube_size^3,...
    mean(cc297_surf_areas)*10000/cube_size^3-prctile(cc297_surf_areas,lower_pct)*10000/cube_size^3,...
    prctile(cc297_surf_areas,upper_pct)*10000/cube_size^3-mean(cc297_surf_areas)*10000/cube_size^3)
hold on
scatter(520,mean(labrador_surf_areas)*10000/cube_size^3)
errorbar(520,mean(labrador_surf_areas)*10000/cube_size^3,...
    mean(labrador_surf_areas)*10000/cube_size^3-prctile(labrador_surf_areas,lower_pct)*10000/cube_size^3,...
    prctile(labrador_surf_areas,upper_pct)*10000/cube_size^3-mean(labrador_surf_areas)*10000/cube_size^3)
hold on
scatter(400,mean(millepora_surf_areas)*10000/cube_size^3)
errorbar(400,mean(millepora_surf_areas)*10000/cube_size^3,...
    mean(millepora_surf_areas)*10000/cube_size^3-prctile(millepora_surf_areas,lower_pct)*10000/cube_size^3,...
    prctile(millepora_surf_areas,upper_pct)*10000/cube_size^3-mean(millepora_surf_areas)*10000/cube_size^3)
hold on
scatter(400,mean(loripes_surf_areas)*10000/cube_size^3)
errorbar(400,mean(loripes_surf_areas)*10000/cube_size^3,...
    mean(loripes_surf_areas)*10000/cube_size^3-prctile(loripes_surf_areas,lower_pct)*10000/cube_size^3,...
    prctile(loripes_surf_areas,upper_pct)*10000/cube_size^3-mean(loripes_surf_areas)*10000/cube_size^3)
hold on
scatter(400,mean(cytherea_surf_areas)*10000/cube_size^3)
errorbar(400,mean(cytherea_surf_areas)*10000/cube_size^3,...
    mean(cytherea_surf_areas)*10000/cube_size^3-prctile(cytherea_surf_areas,lower_pct)*10000/cube_size^3,...
    prctile(cytherea_surf_areas,upper_pct)*10000/cube_size^3-mean(cytherea_surf_areas)*10000/cube_size^3)
hold on
scatter(400,mean(caroliana_surf_areas)*10000/cube_size^3)
errorbar(400,mean(caroliana_surf_areas)*10000/cube_size^3,...
    mean(caroliana_surf_areas)*10000/cube_size^3-prctile(caroliana_surf_areas,lower_pct)*10000/cube_size^3,...
    prctile(caroliana_surf_areas,upper_pct)*10000/cube_size^3-mean(caroliana_surf_areas)*10000/cube_size^3)
hold on
scatter(425,mean(madracis6m_surf_areas)*10000/cube_size^3)
errorbar(425,mean(madracis6m_surf_areas)*10000/cube_size^3,...
    mean(madracis6m_surf_areas)*10000/cube_size^3-prctile(madracis6m_surf_areas,lower_pct)*10000/cube_size^3,...
    prctile(madracis6m_surf_areas,upper_pct)*10000/cube_size^3-mean(madracis6m_surf_areas)*10000/cube_size^3)
hold on
scatter(425,mean(madracis15m_surf_areas)*10000/cube_size^3)
errorbar(425,mean(madracis15m_surf_areas)*10000/cube_size^3,...
    mean(madracis15m_surf_areas)*10000/cube_size^3-prctile(madracis15m_surf_areas,lower_pct)*10000/cube_size^3,...
    prctile(madracis15m_surf_areas,upper_pct)*10000/cube_size^3-mean(madracis15m_surf_areas)*10000/cube_size^3)
hold on
scatter(425,mean(madracis20m_surf_areas)*10000/cube_size^3)
errorbar(425,mean(madracis20m_surf_areas)*10000/cube_size^3,...
    mean(madracis20m_surf_areas)*10000/cube_size^3-prctile(madracis20m_surf_areas,lower_pct)*10000/cube_size^3,...
    prctile(madracis20m_surf_areas,upper_pct)*10000/cube_size^3-mean(madracis20m_surf_areas)*10000/cube_size^3)
xlabel('Time [Ma]')
xlim([350, 650])
ylabel('Surface area / enclosing volume [cm^{-1}]')
set(gca, 'xdir', 'reverse')
%%
figure()
for i = 1:size(c_ord,1)
   scatter(1,i,'filled')
   hold on
end