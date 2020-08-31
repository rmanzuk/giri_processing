function [inner_dense_slices,outer_dense_slices] = densify_slices(inners,outers,densify_factor)
% This function just takes the data from clicking, goes through each
% individual slice, fits a cubic spline to both the inner and outer data,
% and uses that spline to interpolate and add points to the clicked data.
%
% IN
% inners: 1xn_archaeos cell array containing the cell arrays for inner circles 
% created during data collection. To set this variable up, I just call 
% inners = {inner1, inner2, ..., innern}
%
% outers: 1xn_archaeos cell array containing the cell arrays for outer circles 
% created during data collection. To set this variable up, I just call 
% outers = {outer1, outer2, ..., outern}
%
% densify_factor: The multiple degree by which you would like to densify
% the slices. For example, if you would like each slice to have 3x the
% number of points, densify_factor would be 3.
% 
% OUT
% inner_dense_slices: 1xn_archaeos cell array containing the cell arrays
% for inner circle points densified to the desired degree.
%
% outer_dense_slices: 1xn_archaeos cell array containing the cell arrays
% for outer circle points densified to the desired degree.
%
% Ryan A. Manzuk 08/31/2020
    %% begin the function
    % go through each archaeo indi
    inner_dense_slices = {};
        for i = 1:numel(inners)
            % need to go through slice by slice
            for j = 1:numel(inners{i})
                % define this slice in 3d space given the spacing between slices
                this_slice = inners{i}{j};
                % if the current slice isn't empty, get the points and densify
                if ~isempty(this_slice)
                    x = this_slice(:,1);
                    y = this_slice(:,2);
                    % and how many points do we have in the original set?
                    num_points = length(x);
                    % we will analze the spline like the data is circular
                    % make the spline fit
                    curve = cscvn([x';y']);
                    %evaluate to make denser
                    to_evaluate_radians = linspace(0,curve.breaks(end), densify_factor*num_points);
                    output_xy = ppval(curve,to_evaluate_radians)';
                    inner_dense_slices{i}{j} = [output_xy,zeros(size(output_xy,1),1)];
                else
                    inner_dense_slices{i}{j} =[];
                end
            end 
        end

    outer_dense_slices = {};
        for i = 1:numel(outers)
            % need to go through slice by slice
            for j = 1:numel(outers{i})
                % define this slice in 3d space given the spacing between slices
                this_slice = outers{i}{j};
                % if the current slice isn't empty, get the points and densify
                if ~isempty(this_slice)
                    x = this_slice(:,1);
                    y = this_slice(:,2);
                    % and how many points do we have in the original set?
                    num_points = length(x);
                    % we will analze the spline like the data is circular
                    % make the spline fit
                    curve = cscvn([x';y']);
                    %evaluate to make denser
                    to_evaluate_radians = linspace(0,curve.breaks(end), densify_factor*num_points);
                    output_xy = ppval(curve,to_evaluate_radians)';
                    outer_dense_slices{i}{j} = [output_xy,zeros(size(output_xy,1),1)];
                else
                    outer_dense_slices{i}{j} =[];
                end
            end 
        end
end