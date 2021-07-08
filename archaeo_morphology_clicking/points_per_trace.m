% Script to go through and see how many data points exist per branch, per
% slice in all of the data, to optimize densification or de-densification.

%% load pre-processed, cleaned coral data

load('/Users/ryan/Desktop/branch_angle_project/coral_ct_data/millepora_cleaned.mat');
load('/Users/ryan/Desktop/branch_angle_project/coral_ct_data/loripes_cleaned.mat');
load('/Users/ryan/Desktop/branch_angle_project/coral_ct_data/cytherea_cleaned.mat');
load('/Users/ryan/Desktop/branch_angle_project/coral_ct_data/caroliana_cleaned.mat');

%%
% rearrange the coral cell arrays to be 1xn_branches
millepora_reshaped = reshape_coral_cell(millepora_cleaned);
loripes_reshaped = reshape_coral_cell(loripes_cleaned);
cytherea_reshaped = reshape_coral_cell(cytherea_cleaned);
caroliana_reshaped = reshape_coral_cell(caroliana_cleaned);
clear millepora_cleaned
clear loripes_cleaned
clear cytherea_cleaned
clear caroliana_cleaned


% we also need the edges sorted as if going around the circle, not with
% sequential indices
millepora_resorted = sort_outline_points(millepora_reshaped);
loripes_resorted = sort_outline_points(loripes_reshaped);
cytherea_resorted = sort_outline_points(cytherea_reshaped);
caroliana_resorted = sort_outline_points(caroliana_reshaped);
clear millepora_reshaped
clear loripes_reshaped
clear cytherea_reshaped
clear caroliana_reshaped


%%

n_traces = 0;
total_points = 0;
for j = 1:numel(millepora_resorted)
    [points_here,~] = cellfun(@size, millepora_resorted{j}, 'UniformOutput', false);
    total_points = total_points + sum(cell2mat(points_here));
    traces_here = cellfun(@isempty,millepora_resorted{j});
    n_traces = n_traces + sum(~traces_here);
end
pptrace_millepora = total_points/n_traces;

n_traces = 0;
total_points = 0;
for j = 1:numel(loripes_resorted)
    [points_here,~] = cellfun(@size, loripes_resorted{j}, 'UniformOutput', false);
    total_points = total_points + sum(cell2mat(points_here));
    traces_here = cellfun(@isempty,loripes_resorted{j});
    n_traces = n_traces + sum(~traces_here);
end
pptrace_loripes = total_points/n_traces;

n_traces = 0;
total_points = 0;
for j = 1:numel(cytherea_resorted)
    [points_here,~] = cellfun(@size, cytherea_resorted{j}, 'UniformOutput', false);
    total_points = total_points + sum(cell2mat(points_here));
    traces_here = cellfun(@isempty,cytherea_resorted{j});
    n_traces = n_traces + sum(~traces_here);
end
pptrace_cytherea = total_points/n_traces;

n_traces = 0;
total_points = 0;
for j = 1:numel(caroliana_resorted)
    [points_here,~] = cellfun(@size, caroliana_resorted{j}, 'UniformOutput', false);
    total_points = total_points + sum(cell2mat(points_here));
    traces_here = cellfun(@isempty,caroliana_resorted{j});
    n_traces = n_traces + sum(~traces_here);
end
pptrace_caroliana = total_points/n_traces;

n_traces = 0;
total_points = 0;
for j = 1:numel(madracis6m_resorted)
    [points_here,~] = cellfun(@size, madracis6m_resorted{j}, 'UniformOutput', false);
    total_points = total_points + sum(cell2mat(points_here));
    traces_here = cellfun(@isempty,madracis6m_resorted{j});
    n_traces = n_traces + sum(~traces_here);
end
pptrace_madracis6m = total_points/n_traces;

n_traces = 0;
total_points = 0;
for j = 1:numel(madracis15m_resorted)
    [points_here,~] = cellfun(@size, madracis15m_resorted{j}, 'UniformOutput', false);
    total_points = total_points + sum(cell2mat(points_here));
    traces_here = cellfun(@isempty,madracis15m_resorted{j});
    n_traces = n_traces + sum(~traces_here);
end
pptrace_madracis15m = total_points/n_traces;

n_traces = 0;
total_points = 0;
for j = 1:numel(madracis20m_resorted)
    [points_here,~] = cellfun(@size, madracis20m_resorted{j}, 'UniformOutput', false);
    total_points = total_points + sum(cell2mat(points_here));
    traces_here = cellfun(@isempty,madracis20m_resorted{j});
    n_traces = n_traces + sum(~traces_here);
end
pptrace_madracis20m = total_points/n_traces;
%%
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/sm_all_outers.mat');
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/cc297_all_outers.mat');
load('/Users/ryan/Desktop/branch_angle_project/archaeo_clicking_data/labrador_all_outers.mat');

%%
n_traces = 0;
total_points = 0;
for j = 1:numel(sm_outers)
    [points_here,~] = cellfun(@size, sm_outers{j}, 'UniformOutput', false);
    total_points = total_points + sum(cell2mat(points_here));
    traces_here = cellfun(@isempty,sm_outers{j});
    n_traces = n_traces + sum(~traces_here);
end
pptrace_sm = total_points/n_traces;

n_traces = 0;
total_points = 0;
for j = 1:numel(cc297_outers)
    [points_here,~] = cellfun(@size, cc297_outers{j}, 'UniformOutput', false);
    total_points = total_points + sum(cell2mat(points_here));
    traces_here = cellfun(@isempty,cc297_outers{j});
    n_traces = n_traces + sum(~traces_here);
end
pptrace_cc297 = total_points/n_traces;

n_traces = 0;
total_points = 0;
for j = 1:numel(labrador_outers)
    [points_here,~] = cellfun(@size, labrador_outers{j}, 'UniformOutput', false);
    total_points = total_points + sum(cell2mat(points_here));
    traces_here = cellfun(@isempty,labrador_outers{j});
    n_traces = n_traces + sum(~traces_here);
end
pptrace_labrador = total_points/n_traces;