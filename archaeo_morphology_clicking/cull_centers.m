function [cleaned_center_points] = cull_centers(center_points,var_thresh)
    cleaned_center_points = center_points;
    for i = 1:numel(center_points)
        center_spline = cscvn(center_points{i}');
        spline_derivative = fnder(center_spline);
        % set up a set of values to evaluate and do so on the spline and its
        % derivative
        to_eval = linspace(0,center_spline.breaks(end), size(center_points{i},1));
        center_line_eval = ppval(center_spline,to_eval)';
        derivative_eval = ppval(spline_derivative,to_eval)';

        if max(mean(movvar(derivative_eval,3))) > var_thresh
            cleaned_center_points{i} = [NaN,NaN,NaN];
        end
    end
end