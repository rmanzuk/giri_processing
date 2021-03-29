function [outline_dense_slices] = densify_slices(clicked_outlines,densify_factor)
% This function just takes the data from clicking, goes through each
% individual slice, fits a cubic spline to both the inner and outer data,
% and uses that spline to interpolate and add points to the clicked data.
%
% IN
% clicked_outlines: 1xn_branches cell array containing the outlines,
% separated slice-wise, created though automated or clicking data
% collection.
%
% densify_factor: The multiple degree by which you would like to densify
% the slices. For example, if you would like each slice to have 3x the
% number of points, densify_factor would be 3.
% 
% OUT
% outline_dense_slices: 1xn_branches cell array containing slice-wise cell
% arrays with the newly densified outline.
%
%
% Ryan A. Manzuk 08/31/2020
    %% begin the function
    % go through each archaeo indi
    outline_dense_slices = {};
        for i = 1:numel(clicked_outlines)
            % need to go through slice by slice
            for j = 1:numel(clicked_outlines{i})
                % define this slice in 3d space given the spacing between slices
                this_slice = clicked_outlines{i}{j};
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
                    outline_dense_slices{i}{j} = [output_xy,zeros(size(output_xy,1),1)];
                else
                    outline_dense_slices{i}{j} =[];
                end
            end 
        end
end