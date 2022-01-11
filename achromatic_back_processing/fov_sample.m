function [varargout] = fov_sample(varargin)
% this function takes a segmented image and samples it at randomly place
% windows of an input size to test for the impact of field of view on
% modality data.
%
% IN:
% 
% segmented_image: mxn matrix representing the classified image where
% integer values indicate each separate class.
%
% fov_size: 1x2 vector containing the row and column dimensions of the
% desired field of view to be sampled.
%
% n_samples: numer of random fields of view to sample.
%
% optional: 'PointCount' flags that point counting should be done with in
% each window. This input should be followed by a string indicating the
% method ('fixed' for fixed grid, 'stochastic' for randomly placed points,
% and 'both' for both methods). And another input following the method with
% the number of points to be counted. So (....,'PointCount','fixed',100)
% would count 100 points in a fixed, evenly spaced grid pattern in each
% sampled fied of view.
%
% OUT: 
%
% seg_fractions: n_samples x n_classes 2D matrix containing the fraction of
% each class containted withing each randomly sampled field of view of the
% full segmentation.
%
% stoch_point_counts: n_samples x n_classes 2D matrix containing the
% fraction of each class counted with randomly placed points in each sample
% if this was chosen as an input.
%
% fixed_point_counts: n_samples x n_classes 2D matrix containing the
% fraction of each class counted within a fixed grid in each sample
% if this was chosen as an input.
% 
% Ryan A. Manzuk 11/15/2021   
    %% begin the function
    
    % parse the inputs
    segmented_image = varargin{1};
    fov_size = varargin{2};
    n_samples = varargin{3};

    % check if we're point counting in those fields of view
    point_count = false;
    fixed_grid = false;
    stochastic = false;
    if any(strcmp('PointCount',varargin))
        % we're point counting, so indicate that
        point_count = true;
        pc_ind = find(strcmp('PointCount', varargin));
        % gather the number of points (2 inputs after PointCount
        n_points = varargin{pc_ind + 2};
        % check which point counting method
        if strcmp(varargin{pc_ind+1},'fixed')
            fixed_grid = true;
        elseif strcmp(varargin{pc_ind+1},'stochastic')
            stochastic = true;
        elseif strcmp(varargin{pc_ind+1},'both')
            fixed_grid = true;
            stochastic = true;
        else
            error('not a valid point counting method.')
        end
    end

    % if we are point counting, and we're fixed grid, we can just set up
    % those coordinates right now
    if fixed_grid
        % the grid spacing is equal to the square root of the total number
        % of pixels / number of points to be counted.
        grid_spacing = round(sqrt(fov_size(1) * fov_size(2)/n_points));

        % set up the grid. Start with just a count up to the n_points in x
        % direction (columns)
        grid_inds_x = 1:(ceil((fov_size(2)-grid_spacing)/grid_spacing));

        % then just multiply that by a repeated matrix of the values
        % counting the requisite number of times in the y direction (rows).
        grid_inds_x = grid_spacing.*sort(repmat(grid_inds_x,1,(ceil((fov_size(1)-grid_spacing)/grid_spacing))))';
        
        % opposite stuff for y inds
        grid_inds_y = 1:(ceil((fov_size(1)-grid_spacing)/grid_spacing));
        grid_inds_y = grid_spacing.*repmat(grid_inds_y,1,(ceil((fov_size(2)-grid_spacing)/grid_spacing)))';
        
        % and making these into linear inds will be most unseful
        linear_inds = sub2ind(fov_size,grid_inds_y,grid_inds_x);
    end

    % figure out how many classes we're dealing with
    n_classes = numel(unique(segmented_image));

    % pick our random starting points
    fov_corners = round(rand(n_samples,2) .* (size(segmented_image) - fov_size-1));

    % set up empty array for the fractions
    seg_fractions = zeros(n_samples,n_classes);

    % set up arrays for point count fractions if needed
    if stochastic
        stoch_point_fracs = zeros(n_samples,n_classes);
    end
    if fixed_grid
        fixed_point_fracs = zeros(n_samples,n_classes);
    end

    % alright, now just iterate through the number of samplings and get the
    % data
    for i = 1:n_samples

        % extract the fov
        this_bit = segmented_image(fov_corners(i,1):fov_corners(i,1)+fov_size(1)-1,...
        fov_corners(i,2):fov_corners(i,2)+fov_size(2)-1);
        
        % assess the fraction of each class
        for j = 1:n_classes
            pix_present = sum(this_bit == j, 'all');
            seg_fractions(i,j) = pix_present / numel(this_bit);
        end

        % if we're point counting, we have to do that
        if point_count
            % if we want stochastic, do that
            if stochastic
                % set up the sampling points
                point_coords = ceil(rand(n_points,2) .* fov_size);
                % reduce field of view to linear indeces
                linear_inds = sub2ind(fov_size,point_coords(:,1),point_coords(:,2));
                % extract the points
                point_counts = this_bit(linear_inds);
                % make into fractions
                for j = 1:n_classes
                    stoch_point_fracs(i,j) = sum(point_counts == j) / n_points;
                end
            end
            % if we want fixed grid, do that
            if fixed_grid
                % we already have the grid coords, so just skip to sampling
                point_counts = this_bit(linear_inds);
                % make into fractions
                for j = 1:n_classes
                    fixed_point_fracs(i,j) = sum(point_counts == j) / n_points;
                end
            end
        end
    end

    % set up the outputs
    varargout{1} = seg_fractions;
    if stochastic && fixed_grid
        varargout{2} = stoch_point_fracs;
        varargout{3} = fixed_point_fracs;
    elseif stochastic
        varargout{2} = stoch_point_fracs;
    elseif fixed_grid
        varargout{2} = fixed_point_fracs;
    end
end