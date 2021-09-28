% wrapper to read in pbdb data and quick make some diversity curves
%
% R.A. Manzuk 08/09/2021
%% read in the data
raw_data = readtable('/Users/ryan/Desktop/branch_angle_project/pbdb_data/pbdb_data.csv');

taxa_rank = raw_data.accepted_rank;
taxa_name = raw_data.accepted_name;

time_interval = raw_data.early_interval;
min_ages = raw_data.min_ma;
max_ages = raw_data.max_ma;

classes = raw_data.class;

%% filter out the ones that don't have genus level reporting
is_genus = strcmp('genus',taxa_rank) & ~strcmp('Early Cambrian',time_interval) & ~strcmp('Cambrian',time_interval) & ~strcmp('Middle Cambrian',time_interval);

final_taxa = taxa_name(is_genus);
final_time = time_interval(is_genus);
final_min_ages = min_ages(is_genus);
final_max_ages = max_ages(is_genus);

%% distinguish archaeocyathids
is_archaeo = strcmp('Archaeocyatha',classes(is_genus));

archaeo_taxa = final_taxa(is_archaeo);
non_archaeo_taxa = final_taxa(~is_archaeo);

archaeo_times = final_time(is_archaeo);
non_archaeo_times = final_time(~is_archaeo);

archaeo_min_ages = final_min_ages(is_archaeo);
archaeo_max_ages = final_max_ages(is_archaeo);
non_archaeo_min_ages = final_min_ages(~is_archaeo);
non_archaeo_max_ages = final_max_ages(~is_archaeo);

%% assess the number of unique genera in each time bin
non_archaeo_t_bins = unique(non_archaeo_times);
non_archaeo_div = zeros(1,numel(non_archaeo_t_bins));
na_interval_base_ages = zeros(1,numel(non_archaeo_t_bins));
na_interval_top_ages = zeros(1,numel(non_archaeo_t_bins));
for i = 1:numel(non_archaeo_t_bins)
    in_this_bin = strcmp(non_archaeo_t_bins{i},non_archaeo_times);
    unique_genera = unique(non_archaeo_taxa(in_this_bin));
    non_archaeo_div(i) = numel(unique_genera);
    na_interval_base_ages(i) = mean(non_archaeo_min_ages(in_this_bin));
    na_interval_top_ages(i) = mean(non_archaeo_max_ages(in_this_bin));
end

archaeo_t_bins = unique(archaeo_times);
archaeo_div = zeros(1,numel(archaeo_t_bins));
interval_base_ages = zeros(1,numel(archaeo_t_bins));
interval_top_ages = zeros(1,numel(archaeo_t_bins));
for i = 1:numel(archaeo_t_bins)
    in_this_bin = strcmp(archaeo_t_bins{i},archaeo_times);
    unique_genera = unique(archaeo_taxa(in_this_bin));
    archaeo_div(i) = numel(unique_genera);
    interval_base_ages(i) = mean(archaeo_min_ages(in_this_bin));
    interval_top_ages(i) = mean(archaeo_max_ages(in_this_bin));
end

mean_ages = (interval_base_ages + interval_top_ages)/2;
na_mean_ages = (na_interval_base_ages + na_interval_top_ages)/2;

%% re-bin the ages
bin_edges = [492:8:540,650];
na_bin_assignments = discretize(na_mean_ages,bin_edges);
archaeo_bin_assignments = discretize(mean_ages,bin_edges);

rebinned_archaeos = zeros(1,numel(bin_edges)-1);
rebinned_non_archaeos = zeros(1,numel(bin_edges)-1);
for i = 1:numel(bin_edges)-1
    rebinned_archaeos(i) = sum(archaeo_div(archaeo_bin_assignments == i));
    rebinned_non_archaeos(i) = sum(non_archaeo_div(na_bin_assignments == i));
end

bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;

% make it all relative
relative_archaeo_div = rebinned_archaeos/max(rebinned_archaeos);
relative_non_archaeo_div = rebinned_non_archaeos/max(rebinned_non_archaeos);
%% make a figure
figure()
plot(bin_centers,relative_non_archaeo_div)
hold on
plot(bin_centers,relative_archaeo_div)
set(gca,'XDir','reverse')