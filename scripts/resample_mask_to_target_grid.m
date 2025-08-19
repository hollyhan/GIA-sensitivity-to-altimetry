function [mask_resampled] = resample_mask_to_target_grid(mask_data, mask_lat, mask_lon, target_lat, target_lon)
    % resample_mask_to_target_grid - Mask resampling using spatial binning
    %
    % Inputs:
    %   mask_data: binary mask data [lat x lon x time] or [lat x lon] 
    %   mask_lat, mask_lon: source mask coordinates [lat x lon]
    %   target_lat, target_lon: target altimetry coordinates [lat x lon]
    %
    % Outputs:
    %   mask_resampled: resampled mask matching target grid size [lat x lon x time] or [lat x lon]
    %
    % Logic:
    %   For each altimetry pixel at (i,j):
    %   1. Define geographic bounds of that pixel
    %   2. Find ALL glacier mask pixels within those bounds  
    %   3. Count: how many are ice vs. non-ice?
    %   4. Majority vote: >50% ice → output = ice

    % Get dimensions
    if ndims(mask_data) == 3
        [~, ~, ntime] = size(mask_data);
        time_dim = true;
    else
        ntime = 1;
        time_dim = false;
    end
    
    % Get target grid size - use same dimensions as target arrays
    [target_dim1, target_dim2] = size(target_lat);
    
    % Debug: Check target grid dimensions
    fprintf('Target grid dimensions: %dx%d\n', target_dim1, target_dim2);
    
    fprintf('Resampling mask: %dx%d → %dx%d\n', ...
        size(mask_data,1), size(mask_data,2), target_dim1, target_dim2);
    
    % Get target grid bounds and resolution
    target_lat_min = min(target_lat(:));
    target_lat_max = max(target_lat(:));
    target_lon_min = min(target_lon(:));
    target_lon_max = max(target_lon(:));
    
    % Estimate target resolution
    if target_dim1 > 1 && target_dim2 > 1
        target_lat_res = (target_lat_max - target_lat_min) / (target_dim1 - 1);
        target_lon_res = (target_lon_max - target_lon_min) / (target_dim2 - 1);
    else
        % Fallback to interpolation for very small grids
        warning('Target grid too small, falling back to interpolation');
        mask_resampled = resample_mask_to_target_grid(mask_data, mask_lat, mask_lon, target_lat, target_lon);
        return;
    end
    
    % Create target grid vectors
    target_lat_vec = linspace(target_lat_min, target_lat_max, target_dim1);
    target_lon_vec = linspace(target_lon_min, target_lon_max, target_dim2);
    
    % Pre-filter mask data by bounds with buffer
    buffer = max(target_lat_res, target_lon_res) * 2; % 2-pixel buffer
    mask_bounds = mask_lat >= (target_lat_min - buffer) & mask_lat <= (target_lat_max + buffer) & ...
                  mask_lon >= (target_lon_min - buffer) & mask_lon <= (target_lon_max + buffer);
    
    % Get valid mask points
    [valid_rows, valid_cols] = find(mask_bounds);
    if isempty(valid_rows)
        warning('No mask data overlaps with target grid');
        if time_dim
            mask_resampled = false(target_dim1, target_dim2, ntime);
        else
            mask_resampled = false(target_dim1, target_dim2);
        end
        return;
    end
    
    % Extract only valid coordinates (vectorized)
    valid_mask_lat = mask_lat(mask_bounds);
    valid_mask_lon = mask_lon(mask_bounds);
    
    % Find target grid indices for each valid mask point using histcounts
    % Each "cell" corresponds to exactly one pixel in the target altimetry grid
    lat_edges = [target_lat_vec - target_lat_res/2, target_lat_vec(end) + target_lat_res/2];
    lon_edges = [target_lon_vec - target_lon_res/2, target_lon_vec(end) + target_lon_res/2];
    
    % Get bin indices (not values!) for each coordinate
    lat_bins = discretize(valid_mask_lat, lat_edges);
    lon_bins = discretize(valid_mask_lon, lon_edges);
    
    % Ensure both arrays have same length (histcounts can return different sizes)
    min_length = min(length(lat_bins), length(lon_bins));
    lat_bins = lat_bins(1:min_length);
    lon_bins = lon_bins(1:min_length);
    
    % Remove points outside target grid and debug
    fprintf('Bin ranges - lat: %d to %d (target: 1 to %d), lon: %d to %d (target: 1 to %d)\n', ...
        min(lat_bins), max(lat_bins), target_dim1, min(lon_bins), max(lon_bins), target_dim2);
    
    valid_points = lat_bins > 0 & lat_bins <= target_dim1 & lon_bins > 0 & lon_bins <= target_dim2;
    
    % Debug: show how many points are filtered out
    fprintf('Filtering: %d/%d points are within target grid bounds\n', sum(valid_points), length(valid_points));
    
    lat_bins = lat_bins(valid_points);
    lon_bins = lon_bins(valid_points);
    
    % Store the valid points indices for filtering mask data later
    valid_indices_for_mask = false(size(valid_mask_lat));
    valid_indices_for_mask(1:min_length) = valid_points;
    
    if isempty(lat_bins)
        warning('No mask points map to target grid');
        if time_dim
            mask_resampled = false(target_dim1, target_dim2, ntime);
        else
            mask_resampled = false(target_dim1, target_dim2);
        end
        return;
    end
    
    % Final safety check before sub2ind
    fprintf('Final bin ranges - lat: %d to %d, lon: %d to %d\n', ...
        min(lat_bins), max(lat_bins), min(lon_bins), max(lon_bins));
    
    % Convert to linear indices for accumarray
    target_indices = sub2ind([target_dim1, target_dim2], lat_bins, lon_bins);
    
    % Initialize output
    if time_dim
        mask_resampled = false(target_dim1, target_dim2, ntime);
    else
        mask_resampled = false(target_dim1, target_dim2);
    end
    
    % Process each time slice
    for t = 1:ntime
        if time_dim
            current_mask = mask_data(:,:,t);
        else
            current_mask = mask_data;
        end
        
        % Extract valid mask values
        valid_mask_values = current_mask(mask_bounds);
        % Apply the same filtering as coordinates
        valid_mask_values = valid_mask_values(valid_indices_for_mask);
        
        % Use accumarray for ultra-fast majority vote
        % Count ice pixels in each target grid cell
        ice_counts = accumarray(target_indices, double(valid_mask_values), [target_dim1 * target_dim2, 1]);
        total_counts = accumarray(target_indices, 1, [target_dim1 * target_dim2, 1]);
        
        % Majority vote: ice if >50% of pixels in cell are ice
        ice_fraction = ice_counts ./ max(total_counts, 1); % Avoid division by zero
        ice_majority = ice_fraction > 0.5;
        
        % Reshape back to target grid
        result_grid = reshape(ice_majority, target_dim1, target_dim2);
        
        % Fix orientation based on coordinate direction (only for first time slice)
        if t == 1
            % Check if latitude decreases from top to bottom (needs flip)
            lat_decreases = target_lat(1,1) > target_lat(end,1);
            fprintf('Coordinate orientation: lat decreases = %s\n', mat2str(lat_decreases));
        end
        
        % Apply flip only if latitude decreases from top to bottom
        if target_lat(1,1) > target_lat(end,1)
            result_grid = flipud(result_grid);
        end
        
        if time_dim
            mask_resampled(:,:,t) = result_grid;
        else
            mask_resampled = result_grid;
        end
        
        if time_dim
            fprintf('Completed time slice %d/%d\n', t, ntime);
        end
    end
    
    fprintf('Mask resampling completed!\n');
end