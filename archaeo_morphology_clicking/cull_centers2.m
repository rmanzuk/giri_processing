function [cleaned_center_points] = cull_centers2(center_points,length_thresh)
    cleaned_center_points = center_points;
    for i = 1:numel(center_points)
        len = arclength(center_points{i}(:,1),center_points{i}(:,2),center_points{i}(:,3));

        if len < length_thresh
            cleaned_center_points{i} = [NaN,NaN,NaN];
        end
    end
end