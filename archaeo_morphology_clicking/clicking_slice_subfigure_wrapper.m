
% load everything
load('/Users/ryan/Desktop/figures/branch_angle_project/clicking_morphology_method/method_figure_inners.mat');
load('/Users/ryan/Desktop/figures/branch_angle_project/clicking_morphology_method/method_figure_outers.mat');
input_folder = '/Users/ryan/Desktop/figures/branch_angle_project/clicking_morphology_method';
% set up colors
colors = get(gca,'colororder');
% set up image file to read
file_pattern = fullfile(input_folder, '*.tif');
tifs = dir(file_pattern);
fig_image = imread(fullfile(input_folder, tifs(1).name));
% show the image
figure();
imshow(fig_image)
hold on
% iterate through all outlines traced
for j = 1:numel(method_figure_outers)
    
    these_outers = method_figure_outers{j}{1};
    these_inners = method_figure_inners{j}{1};
    % check if outers are empty
    if ~isempty(these_outers)
        % if there are data, define x, y coordinates
        x = these_outers(:,1);
        y = these_outers(:,2);
        % and how many points do we have in the original set?
        num_points = length(x);
        % we will analze the spline like the data is circular
        % make the spline fit
        outer_curve = cscvn([x';y']);
        % put everything on a plot with colors
        scatter(x,y, 20, colors(5,:))
        before = findall(gca);
        fnplt(outer_curve)
        added = setdiff(findall(gca), before);
        set(added, 'Color', colors(6,:))
        set(added, 'LineWidth',1.2)
    end
    % check if inners are empty
    if ~isempty(these_inners)
        % if there are data, define x, y coordinates
        x = these_inners(:,1);
        y = these_inners(:,2);
        % and how many points do we have in the original set?
        num_points = length(x);
        % we will analze the spline like the data is circular
        % make the spline fit
        inner_curve = cscvn([x';y']);
        % put everything on a plot with colors
        scatter(x,y,20, colors(7,:))
        before = findall(gca);
        fnplt(inner_curve)
        added = setdiff(findall(gca), before);
        set(added, 'Color', colors(3,:))
        set(added, 'LineWidth',1.2)
    end
end 