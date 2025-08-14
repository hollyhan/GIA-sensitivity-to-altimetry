function [md_updated, mass_conservation_report] = interpolate_altimetry_to_mesh(h_annual, lat_sphere, long_sphere, years_altimetry, md, X, Y, dhdt_annual, ice_mask)
    % interpolate_altimetry_to_mesh - Interpolates annual ice thickness data onto a global mesh with mass conservation.
    %
    % Usage:
    %   [md_updated, mass_conservation_report] = interpolate_altimetry_to_mesh(h_annual, lat_sphere, long_sphere, years_altimetry, md, X, Y, dhdt_annual, ice_mask)
    %
    % Inputs:
    %   h_annual:         Annual ice thickness data [ny x nx x nt+1] or [nx x ny x nt+1].
    %   lat_sphere:       Latitude coordinates of source data [ny x nx].
    %   long_sphere:      Longitude coordinates of source data [ny x nx].
    %   years_altimetry:  Time vector corresponding to h_annual [nt x 1] or [1 x nt].
    %   md:        Target global mesh model object.
    %   X, Y:             Original altimetry grid coordinates [1 x nx], [1 x ny] - used to verify mass conservation.
    %   dhdt_annual:      Annual elevation change data for mass change verification [ny x nx x nt] - used to verify mass change conservation.
    %
    % Outputs:
    %   md_updated: Updated global mesh model with interpolated thickness in
    %                      md.masstransport.spcthickness.
    %   mass_conservation_report: Structure with mass conservation statistics.

    rhoi = 917.0; % Ice density in kg/m^3

    disp('Interpolating ice altimetry data onto the global mesh...');

    % --- Input Validation ---
    if ~exist('h_annual', 'var') || ~exist('lat_sphere', 'var') || ~exist('long_sphere', 'var') || ~exist('years_altimetry', 'var') || ~exist('md', 'var') || ~exist('X', 'var') || ~exist('Y', 'var') || ~exist('dhdt_annual', 'var')
        error('interpolate_altimetry_to_mesh: Missing required input arguments. h_annual, lat_sphere, long_sphere, years_altimetry, md, X, Y, and dhdt_annual are all required.');
    end
    
    % Mass conservation is always enabled
    enable_mass_conservation = true;
    disp('Mass conservation enabled - interpolating h_annual directly and using dhdt_annual for mass change verification.');
    
    % Check for optional ice mask
    use_ice_mask = false;
    if nargin >= 9 && exist('ice_mask', 'var') && ~isempty(ice_mask)
        use_ice_mask = true;
        disp('Ice mask provided - will restrict interpolation to ice-covered areas.');
    end

    % --- Prepare Data ---
    % Check input sizes and permute h_annual if necessary
    if size(lat_sphere,1) ~= size(long_sphere,1) || size(lat_sphere,2) ~= size(long_sphere,2)
        error('interpolate_altimetry_to_mesh: Source latitude and longitude dimensions must match.');
    end
    if size(h_annual,1) == size(lat_sphere,2) && size(h_annual,2) == size(lat_sphere,1)
        disp('Permuting h_annual dimensions to match lat/long...');
        h_annual = permute(h_annual, [2 1 3]);
    elseif size(h_annual,1) ~= size(lat_sphere,1) || size(h_annual,2) ~= size(lat_sphere,2)
        error('interpolate_altimetry_to_mesh: h_annual spatial dimensions do not match coordinates.');
    end

    [~, ~, nt] = size(h_annual);
    if length(years_altimetry) ~= nt
        error('interpolate_altimetry_to_mesh: Length of years_altimetry must match time dimension of h_annual.');
    end

    % Flatten source coordinates
    lat_source_flat = lat_sphere(:);
    long_source_flat = long_sphere(:);

    % Target mesh info
    n_target_vertices = md.mesh.numberofvertices;
    lat_target = md.mesh.lat;
    long_target = md.mesh.long;
    
    % Initialize the spcthickness field with proper dimensions
    md.masstransport.spcthickness = zeros(n_target_vertices + 1, nt);
    fprintf('Initialized spcthickness field: %d vertices x %d time steps\n', n_target_vertices, nt);
    
    % CRITICAL: Define data coverage bounds to prevent global extrapolation
    % Get actual bounds of altimetry data
    data_lat_min = min(lat_sphere(:));
    data_lat_max = max(lat_sphere(:));
    data_lon_min = min(long_sphere(:));
    data_lon_max = max(long_sphere(:));
    
    fprintf('Altimetry data coverage: lat %.1f to %.1f, lon %.1f to %.1f\n', ...
            data_lat_min, data_lat_max, data_lon_min, data_lon_max);
    
    % Create initial coverage mask using convex hull of actual data points
    % This is more precise than a rectangular bounding box
    valid_data_mask = ~isnan(lat_sphere(:)) & ~isnan(long_sphere(:));
    lat_valid = lat_sphere(valid_data_mask);
    long_valid = long_sphere(valid_data_mask);
    
    % Create convex hull of valid data points
    if length(lat_valid) > 3
        hull_indices = convhull(lat_valid, long_valid);
        hull_lat = lat_valid(hull_indices);
        hull_long = long_valid(hull_indices);
        
        % Check which mesh vertices are inside the convex hull
        data_coverage_mask = inpolygon(lat_target, long_target, hull_lat, hull_long);
        
        fprintf('Using convex hull coverage: %d vertices (%.1f%% of total)\n', ...
                sum(data_coverage_mask), 100*sum(data_coverage_mask)/length(data_coverage_mask));
    else
        % Fallback to rectangular bounds if not enough points
        buffer = 0.1;
        data_lat_min = data_lat_min - buffer;
        data_lat_max = data_lat_max + buffer;
        data_lon_min = data_lon_min - buffer;
        data_lon_max = data_lon_max + buffer;
        
        data_coverage_mask = (lat_target >= data_lat_min) & (lat_target <= data_lat_max) & ...
                             (long_target >= data_lon_min) & (long_target <= data_lon_max);
        
        fprintf('Using rectangular coverage (fallback): %d vertices (%.1f%% of total)\n', ...
                sum(data_coverage_mask), 100*sum(data_coverage_mask)/length(data_coverage_mask));
    end
    
    fprintf('Mesh vertices within data coverage: %d (%.1f%% of total)\n', ...
            sum(data_coverage_mask), 100*sum(data_coverage_mask)/length(data_coverage_mask));

    % Initialize output field in the model
    % h_annual has nt+1 time steps (thickness at start of each year)
    nt_h_annual = size(h_annual, 3);
    nt_dhdt = size(dhdt_annual, 3);
    md.masstransport.spcthickness = zeros(n_target_vertices + 1, nt_h_annual);
    fprintf('Initialized spcthickness field: %d vertices x %d time steps (h_annual)\n', n_target_vertices, nt_h_annual);
    fprintf('Debug: nt_h_annual = %d, nt_dhdt = %d\n', nt_h_annual, nt_dhdt);

    % --- Pre-calculate Mesh Areas (once only) ---
    if enable_mass_conservation
        fprintf('Pre-calculating mesh triangle areas for mass conservation...\n');
        
        % Calculate vertex areas for mass conservation
        vertex_areas = zeros(n_target_vertices, 1);
        
        % Check if mesh has precomputed areas (3D mesh) or needs calculation (2D mesh)
        if isfield(md.mesh, 'area') && ~isempty(md.mesh.area)
            % 3D mesh - use precomputed areas
            triangle_areas = md.mesh.area;
        else
            % 2D mesh - calculate triangle areas manually
            fprintf('  Calculating triangle areas for 2D mesh...\n');
            triangle_areas = zeros(size(md.mesh.elements, 1), 1);
            
            % Debug: Sample some coordinates and calculate area manually (first time only)
            if size(md.mesh.elements, 1) > 0
                elem1 = md.mesh.elements(1, :);
                x1 = md.mesh.x(elem1(1)); y1 = md.mesh.y(elem1(1));
                x2 = md.mesh.x(elem1(2)); y2 = md.mesh.y(elem1(2));
                x3 = md.mesh.x(elem1(3)); y3 = md.mesh.y(elem1(3));
                area1 = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
                
                fprintf('    First triangle vertices: [%d, %d, %d]\n', elem1(1), elem1(2), elem1(3));
                fprintf('    Vertex 1 coords: (%.0f, %.0f)\n', x1, y1);
                fprintf('    Vertex 2 coords: (%.0f, %.0f)\n', x2, y2);
                fprintf('    Vertex 3 coords: (%.0f, %.0f)\n', x3, y3);
                fprintf('    First triangle area: %.0f m² (%.2f km²)\n', area1, area1/1e6);
                
                % Check side lengths
                d12 = sqrt((x2-x1)^2 + (y2-y1)^2);
                d23 = sqrt((x3-x2)^2 + (y3-y2)^2);
                d31 = sqrt((x1-x3)^2 + (y1-y3)^2);
                fprintf('    Triangle side lengths: %.0f, %.0f, %.0f meters\n', d12, d23, d31);

                % Could also use the GetAreas function
                % triangle_areas = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y);
            end
            
            for elem = 1:size(md.mesh.elements, 1)
                % Get triangle vertices
                v1_idx = md.mesh.elements(elem, 1);
                v2_idx = md.mesh.elements(elem, 2);
                v3_idx = md.mesh.elements(elem, 3);
                
                % Use x/y coordinates (projected coordinates in meters, NOT lat/long)
                % This is critical for accurate area calculation in polar stereographic projection
                x1 = md.mesh.x(v1_idx); y1 = md.mesh.y(v1_idx);
                x2 = md.mesh.x(v2_idx); y2 = md.mesh.y(v2_idx);
                x3 = md.mesh.x(v3_idx); y3 = md.mesh.y(v3_idx);
                
                % Calculate triangle area using cross product formula
                triangle_areas(elem) = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
            end
            
            % Debug: Print comprehensive area and size statistics
            min_area = min(triangle_areas);
            max_area = max(triangle_areas);
            mean_area = mean(triangle_areas);
            median_area = median(triangle_areas);
                                        
            fprintf('    === MESH TRIANGLE STATISTICS ===\n');
            fprintf('    Areas     - Min: %.0f m² (%.2f km²), Max: %.0f m² (%.0f km²), Mean: %.0f m² (%.1f km²)\n', ...
                    min_area, min_area/1e6, max_area, max_area/1e6, mean_area, mean_area/1e6);
            fprintf('    Areas     - Median: %.0f m² (%.1f km²)\n', median_area, median_area/1e6);
            % Verify area calculation with alternative method for first triangle
            if size(md.mesh.elements, 1) > 0
                elem1 = md.mesh.elements(1, :);
                x1 = md.mesh.x(elem1(1)); y1 = md.mesh.y(elem1(1));
                x2 = md.mesh.x(elem1(2)); y2 = md.mesh.y(elem1(2));
                x3 = md.mesh.x(elem1(3)); y3 = md.mesh.y(elem1(3));
                
                % Method 1: Cross product (current)
                area1_cross = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
                
                % Method 2: Using vectors
                v1 = [x2-x1, y2-y1];
                v2 = [x3-x1, y3-y1];
                area1_vector = 0.5 * abs(v1(1)*v2(2) - v1(2)*v2(1));
                
                % Method 3: Shoelace formula
                area1_shoelace = 0.5 * abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
                
                fprintf('    AREA CALCULATION VERIFICATION:\n');
                fprintf('      Cross product: %.0f m²\n', area1_cross);
                fprintf('      Vector method: %.0f m²\n', area1_vector);
                fprintf('      Shoelace formula: %.0f m²\n', area1_shoelace);
                fprintf('      Match: %s\n', isequal(round(area1_cross), round(area1_vector), round(area1_shoelace)));
            end
            
            fprintf('    Total mesh area: %.2e m² (%.0f km²)\n', sum(triangle_areas), sum(triangle_areas)/1e6);
            fprintf('    Expected Greenland area: ~2.17e6 km² (for comparison)\n');
            fprintf('    Area ratio (mesh/Greenland): %.2f\n', sum(triangle_areas)/1e6/2.17e6);
            fprintf('    ================================\n');
        end
        
        % Distribute triangle areas to vertices (each triangle contributes 1/3 to each vertex)
        for elem = 1:length(triangle_areas)
            vertices = md.mesh.elements(elem, :);
            for v = 1:3
                if vertices(v) <= n_target_vertices
                    vertex_areas(vertices(v)) = vertex_areas(vertices(v)) + triangle_areas(elem) / 3;
                end
            end
        end
        
        fprintf('Mesh area calculation complete.\n\n');
    end

    % --- Interpolate Each Time Step ---
    fprintf('Interpolating %d time steps...\n', nt);
    
    % Always interpolate h_annual directly to avoid numerical artifacts from small dhdt_annual values
    fprintf('Interpolating h_annual (ice thickness) to mesh...\n');
    data_to_interpolate = h_annual;
    data_type = 'ice thickness';
    

    for t = 1:nt
        % Get current time slice and flatten
        data_slice = data_to_interpolate(:,:,t);
        data_flat = data_slice(:);

        % Identify valid data points for this slice
        valid_mask = ~isnan(lat_source_flat) & ~isnan(long_source_flat) & ~isnan(data_flat);
        
        % Apply ice mask if provided
        if use_ice_mask
            ice_mask_slice = ice_mask(:,:,min(t, size(ice_mask, 3)));
            ice_mask_flat = ice_mask_slice(:);
            valid_mask = valid_mask & ice_mask_flat > 0.5; % Apply ice mask threshold
        end

        if sum(valid_mask) < 3
            warning('interpolate_altimetry_to_mesh: Time step %d: Not enough valid points (%d) for interpolation.', t, sum(valid_mask));
            % Keep NaNs in the output for this time step
            continue;
        end

        % Extract valid source data points
        lat_valid = lat_source_flat(valid_mask);
        long_valid = long_source_flat(valid_mask);
        data_valid = data_flat(valid_mask);

        % Initialize output with zeros (no data outside data coverage)
        data_interp_slice = zeros(n_target_vertices, 1);
        
        % Only interpolate within the data coverage area
        if sum(data_coverage_mask) > 0
            % Create interpolant from valid data with NO extrapolation
            F = scatteredInterpolant(lat_valid, long_valid, data_valid, 'nearest', 'none');
            
            % Interpolate ONLY at mesh points within data coverage
            lat_coverage = lat_target(data_coverage_mask);
            long_coverage = long_target(data_coverage_mask);
            data_coverage = F(lat_coverage, long_coverage);
            
            % Replace NaN values (from failed interpolation) with zeros
            data_coverage(isnan(data_coverage)) = 0;
            
            % Store interpolated values only within coverage area
            data_interp_slice(data_coverage_mask) = data_coverage;
            
            fprintf('  Interpolated at %d vertices (%.1f%% of mesh), %.1f%% have %s > 0\n', ...
                    length(data_coverage), 100*length(data_coverage)/n_target_vertices, ...
                    100*sum(abs(data_coverage) > 1e-6)/length(data_coverage), data_type);
        else
            fprintf('  WARNING: No mesh vertices within data coverage area!\n');
        end

        % Store the interpolated values in the model structure
        md.masstransport.spcthickness(1:n_target_vertices, t) = data_interp_slice;
    end
    
    % --- Mass Conservation Check and Correction ---
    if enable_mass_conservation
        fprintf('\n=== Mass Conservation Analysis ===\n');
        
        % Calculate original mass from altimetry grid
        dx = abs(diff(X(1:2))); % Grid spacing in x-direction (meters)
        dy = abs(diff(Y(1:2))); % Grid spacing in y-direction (meters)
        cell_area_original = dx * dy; % Area of each cell in m^2
        
        % Debug: Print coordinate ranges and grid info
        fprintf('DEBUG - Original grid info:\n');
        fprintf('  X range: %.0f to %.0f meters, dx = %.0f m\n', min(X), max(X), dx);
        fprintf('  Y range: %.0f to %.0f meters, dy = %.0f m\n', min(Y), max(Y), dy);
        fprintf('  Cell area: %.0f m²\n', cell_area_original);
        
        % Debug: Print mesh coordinate info  
        fprintf('DEBUG - Mesh info:\n');
        fprintf('  X range: %.0f to %.0f\n', min(md.mesh.x), max(md.mesh.x));
        fprintf('  Y range: %.0f to %.0f\n', min(md.mesh.y), max(md.mesh.y));
        fprintf('  Number of vertices: %d\n', n_target_vertices);
        fprintf('  Number of triangles: %d\n', size(md.mesh.elements, 1));
        
        % Initialize mass conservation report (only for first nt time steps where mass conservation applies)
        mass_conservation_report.original_mass = zeros(nt, 1);
        mass_conservation_report.interpolated_mass = zeros(nt, 1);
        mass_conservation_report.mass_error = zeros(nt, 1);
        mass_conservation_report.correction_factor = zeros(nt, 1);
        
        for t = 1:nt_h_annual
            fprintf('Debug: Processing time step t = %d (nt_h_annual = %d, nt = %d)\n', t, nt_h_annual, nt);
            % For mass conservation, we only have dhdt_annual data for the first nt time steps
            if t <= nt
                % Calculate mass changes from original dhdt_annual data
                % dhdt_annual(t) represents change from year t to t+1
                mass_change_original = dhdt_annual(:,:,t) * cell_area_original * rhoi * 1e-12; % Convert to Gt
                original_total_mass_change = sum(mass_change_original(:), 'omitnan');
                mass_type = 'annual mass change';
            else
                % For the last time step (t = nt+1), no mass conservation needed
                original_total_mass_change = 0;
                mass_type = 'final thickness';
            end
            
            % Interpolated mass calculation (from whatever we interpolated)
            data_interp_t = md.masstransport.spcthickness(1:n_target_vertices, t);
            
            % Calculate mass on mesh using triangle areas
            % For each vertex, we need to distribute its area among connected triangles
            vertex_areas = zeros(n_target_vertices, 1);
            
            % Check if mesh has precomputed areas (3D mesh) or needs calculation (2D mesh)
            if isfield(md.mesh, 'area') && ~isempty(md.mesh.area)
                % 3D mesh - use precomputed areas
                triangle_areas = md.mesh.area;
            else
                % 2D mesh - calculate triangle areas manually
                fprintf('  Calculating triangle areas for 2D mesh...\n');
                triangle_areas = zeros(size(md.mesh.elements, 1), 1);
                
                % Debug: Sample some coordinates and calculate area manually
                if size(md.mesh.elements, 1) > 0
                    elem1 = md.mesh.elements(1, :);
                    x1 = md.mesh.x(elem1(1)); y1 = md.mesh.y(elem1(1));
                    x2 = md.mesh.x(elem1(2)); y2 = md.mesh.y(elem1(2));
                    x3 = md.mesh.x(elem1(3)); y3 = md.mesh.y(elem1(3));
                    area1 = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
                    
                    fprintf('    First triangle vertices: [%d, %d, %d]\n', elem1(1), elem1(2), elem1(3));
                    fprintf('    Vertex 1 coords: (%.0f, %.0f)\n', x1, y1);
                    fprintf('    Vertex 2 coords: (%.0f, %.0f)\n', x2, y2);
                    fprintf('    Vertex 3 coords: (%.0f, %.0f)\n', x3, y3);
                    fprintf('    First triangle area: %.0f m² (%.2f km²)\n', area1, area1/1e6);
                    
                    % Check side lengths
                    d12 = sqrt((x2-x1)^2 + (y2-y1)^2);
                    d23 = sqrt((x3-x2)^2 + (y3-y2)^2);
                    d31 = sqrt((x1-x3)^2 + (y1-y3)^2);
                    fprintf('    Triangle side lengths: %.0f, %.0f, %.0f meters\n', d12, d23, d31);

                    % Could also use the GetAreas function
                    % triangle_areas = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y);
                end
                
                for elem = 1:size(md.mesh.elements, 1)
                    % Get triangle vertices
                    v1_idx = md.mesh.elements(elem, 1);
                    v2_idx = md.mesh.elements(elem, 2);
                    v3_idx = md.mesh.elements(elem, 3);
                    
                    % Use x/y coordinates (projected coordinates in meters, NOT lat/long)
                    % This is critical for accurate area calculation in polar stereographic projection
                    x1 = md.mesh.x(v1_idx); y1 = md.mesh.y(v1_idx);
                    x2 = md.mesh.x(v2_idx); y2 = md.mesh.y(v2_idx);
                    x3 = md.mesh.x(v3_idx); y3 = md.mesh.y(v3_idx);
                    
                    % Calculate triangle area using cross product formula
                    triangle_areas(elem) = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
                end
                
                % Debug: Print comprehensive area and size statistics
                min_area = min(triangle_areas);
                max_area = max(triangle_areas);
                mean_area = mean(triangle_areas);
                median_area = median(triangle_areas);
                                    
                fprintf('    === MESH TRIANGLE STATISTICS ===\n');
                fprintf('    Areas     - Min: %.0f m² (%.2f km²), Max: %.0f m² (%.0f km²), Mean: %.0f m² (%.1f km²)\n', ...
                        min_area, min_area/1e6, max_area, max_area/1e6, mean_area, mean_area/1e6);
                fprintf('    Areas     - Median: %.0f m² (%.1f km²)\n', median_area, median_area/1e6);
                % Verify area calculation with alternative method for first triangle
                if size(md.mesh.elements, 1) > 0
                    elem1 = md.mesh.elements(1, :);
                    x1 = md.mesh.x(elem1(1)); y1 = md.mesh.y(elem1(1));
                    x2 = md.mesh.x(elem1(2)); y2 = md.mesh.y(elem1(2));
                    x3 = md.mesh.x(elem1(3)); y3 = md.mesh.y(elem1(3));
                    
                    % Method 1: Cross product (current)
                    area1_cross = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
                    
                    % Method 2: Using vectors
                    v1 = [x2-x1, y2-y1];
                    v2 = [x3-x1, y3-y1];
                    area1_vector = 0.5 * abs(v1(1)*v2(2) - v1(2)*v2(1));
                    
                    % Method 3: Shoelace formula
                    area1_shoelace = 0.5 * abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
                    
                    fprintf('    AREA CALCULATION VERIFICATION:\n');
                    fprintf('      Cross product: %.0f m²\n', area1_cross);
                    fprintf('      Vector method: %.0f m²\n', area1_vector);
                    fprintf('      Shoelace formula: %.0f m²\n', area1_shoelace);
                    fprintf('      Match: %s\n', isequal(round(area1_cross), round(area1_vector), round(area1_shoelace)));
                end
                
                fprintf('    Total mesh area: %.2e m² (%.0f km²)\n', sum(triangle_areas), sum(triangle_areas)/1e6);
                fprintf('    Expected Greenland area: ~2.17e6 km² (for comparison)\n');
                fprintf('    Area ratio (mesh/Greenland): %.2f\n', sum(triangle_areas)/1e6/2.17e6);
                fprintf('    ================================\n');
            end
            
            % Distribute triangle areas to vertices (each triangle contributes 1/3 to each vertex)
            for elem = 1:length(triangle_areas)
                vertices = md.mesh.elements(elem, :);
                for v = 1:3
                    if vertices(v) <= n_target_vertices
                        vertex_areas(vertices(v)) = vertex_areas(vertices(v)) + triangle_areas(elem) / 3;
                    end
                end
            end
            
            % Calculate mass changes from interpolated h_annual data
            if t < nt_h_annual
                % h_annual has nt+1 time steps, dhdt_annual has nt time steps
                % For mass change from t to t+1, we need h(t) and h(t+1)
                h_interp_t = data_interp_t;  % This is h(t) - thickness at start of year t
                h_interp_t_plus_1 = md.masstransport.spcthickness(1:n_target_vertices, t+1);  % This is h(t+1) - thickness at start of year t+1
                
                % Calculate mass change: (h(t+1) - h(t)) * area * density
                mass_change_interp_vertices = (h_interp_t_plus_1 - h_interp_t) .* vertex_areas * rhoi * 1e-12; % Convert to Gt
                interpolated_total_mass_change = sum(mass_change_interp_vertices, 'omitnan');
            else
                % For the last time step (t = nt_h_annual), no mass change to calculate
                interpolated_total_mass_change = 0;
            end
            
            % Debug: Check interpolated data distribution
            non_zero_vertices = sum(abs(data_interp_t) > 1e-6);
            max_data = max(abs(data_interp_t));
            mean_data_nonzero = mean(abs(data_interp_t(abs(data_interp_t) > 1e-6)));
            fprintf('    %s stats: %d vertices with data (%.1f%%), max=%.1f, mean=%.1f\n', ...
                    data_type, non_zero_vertices, 100*non_zero_vertices/n_target_vertices, max_data, mean_data_nonzero);
            
            % Calculate mass change error and correction factor (only for first nt time steps)
            if t <= nt
                mass_error = abs(interpolated_total_mass_change - original_total_mass_change);
                if interpolated_total_mass_change ~= 0
                    correction_factor = original_total_mass_change / interpolated_total_mass_change;
                else
                    correction_factor = 1.0;
                end
                
                % Debug: Print mass change calculation details
                fprintf('  DEBUG - Mass change details for time step %d (%s):\n', t, mass_type);
                fprintf('    Original mass change: %.2f Gt\n', original_total_mass_change);
                fprintf('    Interpolated mass change: %.2f Gt\n', interpolated_total_mass_change);
                fprintf('    Error: %.2f Gt (%.1f%%)\n', mass_error, 100*mass_error/abs(original_total_mass_change));
                fprintf('    Correction factor: %.4f\n', correction_factor);
                
                % Apply mass conservation correction
                md.masstransport.spcthickness(1:n_target_vertices, t) = ...
                    md.masstransport.spcthickness(1:n_target_vertices, t) * correction_factor;
            else
                % For the last time step, no mass conservation needed
                correction_factor = 1.0;
                mass_error = 0;
                fprintf('  Time step %d: Final thickness (no mass conservation)\n', t);
            end
            
            % Calculate corrected mass change to verify perfect conservation
            if t < nt_h_annual
                data_corrected_t = md.masstransport.spcthickness(1:n_target_vertices, t);
                h_corrected_t = data_corrected_t;
                h_corrected_t_plus_1 = md.masstransport.spcthickness(1:n_target_vertices, t+1);
                mass_change_corrected_vertices = (h_corrected_t_plus_1 - h_corrected_t) .* vertex_areas * rhoi * 1e-12; % Convert to Gt
                corrected_total_mass_change = sum(mass_change_corrected_vertices, 'omitnan');
            else
                % For the last time step (t = nt_h_annual), no mass change to calculate
                corrected_total_mass_change = 0;
            end
            
            % Store statistics (only for first nt time steps where mass conservation applies)
            if t <= nt
                mass_conservation_report.original_mass(t) = original_total_mass_change;
                mass_conservation_report.interpolated_mass(t) = interpolated_total_mass_change;
                mass_conservation_report.corrected_mass(t) = corrected_total_mass_change;
                mass_conservation_report.mass_error(t) = mass_error;
                mass_conservation_report.correction_factor(t) = correction_factor;
                mass_conservation_report.final_error(t) = abs(corrected_total_mass_change - original_total_mass_change);
            end
            
            fprintf('Time step %d (%s): Original=%.2f Gt, Interpolated=%.2f Gt, Corrected=%.2f Gt\n', ...
                    t, mass_type, original_total_mass_change, interpolated_total_mass_change, corrected_total_mass_change);
            fprintf('              Error before=%.2f Gt (%.1f%%), Error after=%.2e Gt (%.2e%%), Correction=%.4f\n', ...
                    mass_error, 100*mass_error/abs(original_total_mass_change), ...
                    abs(corrected_total_mass_change - original_total_mass_change), ...
                    100*abs(corrected_total_mass_change - original_total_mass_change)/abs(original_total_mass_change), ...
                    correction_factor);
        end
        
        % Summary statistics
        avg_error_before_percent = mean(100 * mass_conservation_report.mass_error ./ abs(mass_conservation_report.original_mass));
        max_error_before_percent = max(100 * mass_conservation_report.mass_error ./ abs(mass_conservation_report.original_mass));
        avg_error_after_percent = mean(100 * mass_conservation_report.final_error ./ abs(mass_conservation_report.original_mass));
        max_error_after_percent = max(100 * mass_conservation_report.final_error ./ abs(mass_conservation_report.original_mass));
        avg_correction = mean(mass_conservation_report.correction_factor);
        
        mass_conservation_report.avg_error_before_percent = avg_error_before_percent;
        mass_conservation_report.max_error_before_percent = max_error_before_percent;
        mass_conservation_report.avg_error_after_percent = avg_error_after_percent;
        mass_conservation_report.max_error_after_percent = max_error_after_percent;
        mass_conservation_report.avg_correction_factor = avg_correction;
        
        % Conservation type for summary
        conservation_type = 'Annual Mass Change Conservation';
        
        fprintf('\n--- %s Summary ---\n', conservation_type);
        fprintf('BEFORE correction - Avg error: %.2f%%, Max error: %.2f%%\n', avg_error_before_percent, max_error_before_percent);
        fprintf('AFTER correction  - Avg error: %.2e%%, Max error: %.2e%%\n', avg_error_after_percent, max_error_after_percent);
        fprintf('Average correction factor: %.4f\n', avg_correction);
        fprintf('=================================\n\n');
    end

    % Only reconstruct if we interpolated elevation changes (not needed for direct h_annual interpolation)
    if strcmp(data_type, 'elevation change')
        fprintf('\n=== Reconstructing Ice Thickness Time Series ===\n');
        fprintf('Converting interpolated dhdt_annual to h_annual time series...\n');
        
        % Create new thickness field by working backwards from assumed final thickness
        h_reconstructed = zeros(n_target_vertices, nt);
        
        % Set final thickness to a reasonable estimate (using mean from original data if available)
        if ~isempty(h_annual)
            % Use mean final thickness from original data where ice exists
            original_final = h_annual(:,:,end);
            mean_final_thickness = mean(original_final(original_final > 0), 'omitnan');
            if isnan(mean_final_thickness) || mean_final_thickness <= 0
                mean_final_thickness = 1500; % Default fallback
            end
        else
            mean_final_thickness = 1500; % Default mean Greenland ice thickness
        end
        
        % Apply this mean thickness only where we have interpolated data
        final_dhdt = md.masstransport.spcthickness(1:n_target_vertices, end);
        has_ice_mask = abs(final_dhdt) > 1e-6; % Where we have actual elevation change data
        
        h_reconstructed(:, end) = 0; % Initialize to zero everywhere (no ice outside coverage)
        h_reconstructed(has_ice_mask, end) = mean_final_thickness; % Set thickness only where ice exists
        
        fprintf('Using mean final thickness of %.1f m for %d vertices with ice data\n', ...
                mean_final_thickness, sum(has_ice_mask));
        
        % Work backwards through time using dhdt to reconstruct thickness
        for t = nt-1:-1:1
            dhdt_t = md.masstransport.spcthickness(1:n_target_vertices, t);
            h_reconstructed(:, t) = h_reconstructed(:, t+1) - dhdt_t; % Subtract dhdt to go backwards in time
            
            % Ensure thickness doesn't go negative where we expect ice
            h_reconstructed(h_reconstructed(:, t) < 0 & has_ice_mask, t) = 0;
        end
        
        % Replace the dhdt data with reconstructed thickness
        md.masstransport.spcthickness(1:n_target_vertices, :) = h_reconstructed;
        
        fprintf('Ice thickness time series reconstructed from elevation changes.\n');
        fprintf('Final thickness range: %.1f to %.1f m (mean: %.1f m for ice-covered vertices)\n', ...
                min(h_reconstructed(has_ice_mask, end)), max(h_reconstructed(has_ice_mask, end)), ...
                mean(h_reconstructed(has_ice_mask, end)));
        
        disp('Interpolation complete - reconstructed h_annual (ice thickness) stored in md.masstransport.spcthickness.');
    else
        disp('Interpolation complete - h_annual (ice thickness) stored in md.masstransport.spcthickness.');
    end

    % Add the time array in the last row of spcthickness
    md.masstransport.spcthickness(end, :) = years_altimetry(:)'; % Ensure time is a row vector

    % Return the updated model
    md_updated = md;
    disp('====================================');
end