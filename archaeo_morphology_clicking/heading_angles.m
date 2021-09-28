function [angles] = heading_angles(outer_center_stats)
% This function calculates the set of angles between all branches and the
% mean heading of the colony.
%
% IN
% outer_center_stats: Structure containing 1xn_archaeo cell arrays where each 
% cell contains a given measure of the outer tube for the given archaeo.
% (output of the center_line_analysis.m function)
%
% OUT
% angles: list of angles between each branch and the mean heading of the
% colony.
%
%
% R. A. Manzuk 09/21/2021
    %% begin the function
    % first, we need to extract all branch derivatives and find their mean
    mean_branch_heading = [];
    for i = 1:numel(outer_center_stats.derivative)
        mean_branch_heading(i,:) = mean(outer_center_stats.derivative{i},'omitnan');
    end
    mean_deriv = mean(mean_branch_heading,'omitnan');
    angles = [];
    for i = 1:size(mean_branch_heading,1)
        this_angle = rad2deg(atan2(norm(cross(mean_deriv,mean_branch_heading(i,:))),dot(mean_deriv,mean_branch_heading(i,:))));
        if ~isnan(this_angle)
            angles = [angles, this_angle];
        end
    end
end