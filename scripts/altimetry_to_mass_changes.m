function [mass_changes, grid_areas, cumulative_mass] = altimetry_to_mass_changes( ...
    dhdt_data, lat_grid, lon_grid, years, start_year, end_year, options)
% Convert altimetry elevation changes to mass changes and compute cumulative totals.
%
% Inputs:
%   dhdt_data - elevation change rates [m/year] (n_lat x n_lon x n_time) or (n_points x n_time)
%   lat_grid, lon_grid - coordinates [degrees]
%   years - vector of time stamps corresponding to 3rd dimension of dhdt_data
%   start_year, end_year - define time window for cumulative integration
%   options - structure with parameters:
%       .rho_ice - ice density [kg/m³] (default: 917)
%       .mask - optional logical mask
%
% Outputs:
%   mass_changes  - instantaneous mass change rates [kg/yr]
%   grid_areas    - grid cell areas [m²]
%   cumulative_mass - cumulative mass change [kg] relative to start_year

    % Parse options
    if nargin < 4 || isempty(options)
        options = struct();
    end
    
    if ~isfield(options, 'rho_ice')
        options.rho_ice = 917; % kg/m³ (ice density)
    end
    
    % Determine data dimensions
    dhdt_size = size(dhdt_data);
    if length(dhdt_size) == 3
        % 3D grid data (lat x lon x time) - e.g., MEaSUREs 1421x783x31
        [n_lat, n_lon, n_time] = size(dhdt_data);
        is_gridded = true;
        
        fprintf('Processing 3D gridded data: %dx%dx%d\n', n_lat, n_lon, n_time);
        
        % Create coordinate grids if vectors provided
        if isvector(lat_grid) && isvector(lon_grid)
            [lon_grid, lat_grid] = meshgrid(lon_grid, lat_grid);
        end
        
        % Validate grid dimensions match data
        if ~isequal(size(lat_grid), [n_lat, n_lon]) || ~isequal(size(lon_grid), [n_lat, n_lon])
            error('Coordinate grid dimensions [%dx%d] do not match data dimensions [%dx%d]', ...
                size(lat_grid), n_lat, n_lon);
        end
        
    elseif length(dhdt_size) == 2
        % 2D scattered data (points x time)
        [n_points, n_time] = size(dhdt_data);
        is_gridded = false;
        
        if ~isvector(lat_grid) || ~isvector(lon_grid)
            error('For 2D dhdt_data, lat_grid and lon_grid must be vectors');
        end
        
        if length(lat_grid) ~= n_points || length(lon_grid) ~= n_points
            error('lat_grid and lon_grid length must match dhdt_data first dimension');
        end
        
    else
        error('dhdt_data must be 2D (points x time) or 3D (lat x lon x time)');
    end

    % Calculate grid areas
    % -------------------
    disp('Calculating grid cell areas...');
    
    if is_gridded
        % For gridded data, calculate area of each grid cell
        grid_areas = calculate_grid_areas_3d(lat_grid, lon_grid);
        
        % Reshape data and coordinates to vectors
        lat_vec = lat_grid(:);
        lon_vec = lon_grid(:);
        dhdt_vec = reshape(dhdt_data, [], n_time); % (n_lat*n_lon) x n_time
        grid_areas = grid_areas(:);
        
        % Apply mask if provided
        if isfield(options, 'mask') && ~isempty(options.mask)
            mask_vec = options.mask(:);
            valid_idx = mask_vec > 0;
            
            lat_vec = lat_vec(valid_idx);
            lon_vec = lon_vec(valid_idx);
            dhdt_vec = dhdt_vec(valid_idx, :);
            grid_areas = grid_areas(valid_idx);
        end
        
    else
        % For scattered data, calculate exact areas from coordinate differences
        fprintf('Calculating exact grid areas from coordinate data...\n');
        grid_areas = calculate_exact_point_areas(lat_grid, lon_grid);
        dhdt_vec = dhdt_data;
    end

    % Convert elevation changes to mass changes
    % ----------------------------------------
    disp('Converting elevation changes to mass changes...');
    
    % mass_change = dhdt * rho_ice * area
    % Expand grid_areas to match time dimension
    areas_expanded = repmat(grid_areas, 1, n_time);
    
    % Calculate mass changes
    mass_changes = dhdt_vec .* options.rho_ice .* areas_expanded / 1e12; % Gt
    
    fprintf('Converted %d data points with %d time steps\n', size(mass_changes, 1), n_time);
    fprintf('Mass change range: %.2e to %.2e Gt\n', min(mass_changes(:)), max(mass_changes(:)));

     % --- Compute cumulative mass change ---
    if nargin >= 6 && ~isempty(start_year) && ~isempty(end_year)
        [~, i_start] = min(abs(years - start_year));
        [~, i_end] = min(abs(years - end_year));
        if i_start >= i_end
            error('end_year must be greater than start_year.');
        end
        size(dhdt_vec)
        % Integrate (sum) dh/dt * rho * area between start and end
        dh_cumulative = cumsum(dhdt_vec(:,1:end-1), 2, 'omitnan');
        size(dh_cumulative)
        cumulative_mass = dh_cumulative(:,i_end-1) * options.rho_ice .* grid_areas / 1e12; % Gt;
    else
        cumulative_mass = [];
    end

    fprintf('Processed %d cells, %d timesteps (%g–%g)\n', size(mass_changes,1), n_time, years(1), years(end));
end


function areas = calculate_grid_areas_3d(lat_grid, lon_grid)
    % Calculate exact area of each grid cell using actual coordinate spacing
    
    [n_lat, n_lon] = size(lat_grid);
    areas = zeros(n_lat, n_lon);
    
    % Earth radius
    R = 6371000; % meters
    
    fprintf('Calculating exact grid cell areas from coordinate data...\n');
    
    for i = 1:n_lat
        for j = 1:n_lon
            % Get actual grid cell boundaries from coordinate data
            
            % Latitude boundaries
            if i == 1
                if n_lat > 1
                    lat_south = lat_grid(i, j) - (lat_grid(2, j) - lat_grid(1, j))/2;
                else
                    lat_south = lat_grid(i, j) - 0.025; % Default small spacing
                end
            else
                lat_south = (lat_grid(i-1, j) + lat_grid(i, j))/2;
            end
            
            if i == n_lat
                if n_lat > 1
                    lat_north = lat_grid(i, j) + (lat_grid(i, j) - lat_grid(i-1, j))/2;
                else
                    lat_north = lat_grid(i, j) + 0.025; % Default small spacing
                end
            else
                lat_north = (lat_grid(i, j) + lat_grid(i+1, j))/2;
            end
            
            % Longitude boundaries
            if j == 1
                if n_lon > 1
                    lon_west = lon_grid(i, j) - (lon_grid(i, 2) - lon_grid(i, 1))/2;
                else
                    lon_west = lon_grid(i, j) - 0.025; % Default small spacing
                end
            else
                lon_west = (lon_grid(i, j-1) + lon_grid(i, j))/2;
            end
            
            if j == n_lon
                if n_lon > 1
                    lon_east = lon_grid(i, j) + (lon_grid(i, j) - lon_grid(i, j-1))/2;
                else
                    lon_east = lon_grid(i, j) + 0.025; % Default small spacing
                end
            else
                lon_east = (lon_grid(i, j) + lon_grid(i, j+1))/2;
            end
            
            % Calculate area using spherical geometry
            lat_south_rad = deg2rad(lat_south);
            lat_north_rad = deg2rad(lat_north);
            lon_west_rad = deg2rad(lon_west);
            lon_east_rad = deg2rad(lon_east);
            
            areas(i, j) = R^2 * (lon_east_rad - lon_west_rad) * ...
                         (sin(lat_north_rad) - sin(lat_south_rad));
        end
    end
    
    fprintf('Grid cell areas: %.0f to %.0f m² (%.1f to %.1f km²)\n', ...
        min(areas(:)), max(areas(:)), min(areas(:))/1e6, max(areas(:))/1e6);
end

function areas = calculate_exact_point_areas(lat_points, lon_points)
    % Calculate exact area for each point using Voronoi-like approach with actual coordinate spacing
    
    n_points = length(lat_points);
    areas = zeros(n_points, 1);
    
    fprintf('Calculating exact point areas from coordinate spacing...\n');
    
    % Earth radius
    R = 6371000; % meters
    
    % For large datasets, use efficient nearest neighbor approach
    if n_points < 10000
        % Exact Voronoi-like calculation for smaller datasets
        for i = 1:n_points
            % Find nearest neighbors in each direction
            lat_i = lat_points(i);
            lon_i = lon_points(i);
            
            % Calculate distances to all other points
            lat_diff = lat_points - lat_i;
            lon_diff = lon_points - lon_i;
            
            % Convert to approximate meters for distance calculation
            mean_lat = lat_i;
            m_per_deg_lat = 111320; % meters per degree of latitude on Earth
            m_per_deg_lon = 111320 * cos(deg2rad(mean_lat));
            
            distances = sqrt((lat_diff * m_per_deg_lat).^2 + (lon_diff * m_per_deg_lon).^2);
            distances(i) = Inf; % Exclude self
            
            % Find nearest neighbors
            [sorted_dist, sorted_idx] = sort(distances);
            
            % Use first few neighbors to estimate grid spacing
            n_neighbors = min(8, n_points-1);
            neighbor_distances = sorted_dist(1:n_neighbors);
            
            % Estimate cell size as half the mean nearest neighbor distance
            mean_spacing_m = mean(neighbor_distances) / 2;
            areas(i) = mean_spacing_m^2;
        end
    else
        % Fast approximation for large datasets
        % Sample a subset to estimate typical grid spacing
        sample_size = min(1000, n_points);
        sample_idx = round(linspace(1, n_points, sample_size));
        
        spacings = zeros(sample_size, 1);
        for k = 1:sample_size
            i = sample_idx(k);
            lat_i = lat_points(i);
            lon_i = lon_points(i);
            
            % Calculate distances to all other points
            lat_diff = lat_points - lat_i;
            lon_diff = lon_points - lon_i;
            
            % Convert to meters
            m_per_deg_lat = 111320;
            m_per_deg_lon = 111320 * cos(deg2rad(lat_i));
            
            distances = sqrt((lat_diff * m_per_deg_lat).^2 + (lon_diff * m_per_deg_lon).^2);
            distances(i) = Inf;
            
            % Find nearest neighbor
            spacings(k) = min(distances);
        end
        
        % Use median spacing for all points
        typical_spacing = median(spacings) / 2;
        areas(:) = typical_spacing^2;
        
        fprintf('Using median grid spacing: %.0f m\n', typical_spacing*2);
    end
    
    fprintf('Point areas: %.0f to %.0f m² (%.1f to %.1f km²)\n', ...
        min(areas), max(areas), min(areas)/1e6, max(areas)/1e6);
end

function areas = estimate_point_areas(lat_points, lon_points)
    % Legacy function - kept for compatibility
    areas = calculate_exact_point_areas(lat_points, lon_points);
end