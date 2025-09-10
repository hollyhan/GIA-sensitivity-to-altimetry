function [md_updated, mass_conservation_report] = interpolate_altimetry_to_mesh(h_annual, lat_sphere, long_sphere, years_dhdt, md, X, Y, dhdt_annual)
    % interpolate_altimetry_to_mesh - Interpolates annual ice thickness data onto a global mesh with mass conservation.
    %
    % Usage:
    %   [md_updated, mass_conservation_report] = interpolate_altimetry_to_mesh(h_annual, lat_sphere, long_sphere, years_dhdt, md, X, Y, dhdt_annual)
    %
    % Inputs:
    %   h_annual:         Annual ice thickness data [ny x nx x nt] or [nx x ny x nt+1].
    %   lat_sphere:       Latitude coordinates of source data [ny x nx].
    %   long_sphere:      Longitude coordinates of source data [ny x nx].
    %   years_dhdt:       Time vector corresponding to dhdt_annual [nt x 1] or [1 x nt].
    %   md:        Target global mesh model object.
    %   X, Y:             Original altimetry grid coordinates [1 x nx], [1 x ny] - REQUIRED for mass conservation.
    %   dhdt_annual:      Annual elevation change data for mass conservation [ny x nx x nt].
    %
    % Outputs:
    %   md_updated: Updated mesh model with interpolated thickness in
    %                      md.masstransport.spcthickness.
    %   mass_conservation_report: Structure with mass conservation statistics.
    %
    % Notes: Recommended to use dhdt_annual for mass conservation.

    rhoi = 917.0; % Ice density in kg/m^3

    disp('Interpolating ice altimetry data onto the global mesh...');

    % --- Input Validation ---
    if ~exist('h_annual', 'var') || ~exist('lat_sphere', 'var') || ~exist('long_sphere', 'var') || ~exist('years_dhdt', 'var') || ~exist('md', 'var') || ~exist('X', 'var') || ~exist('Y', 'var')
        error('interpolate_altimetry_to_mesh: Missing required input arguments. X and Y are required for mass conservation.');
    end
    
    % Mass conservation is always enabled (X, Y required)
    enable_mass_conservation = true;
    if nargin >= 8 && exist('dhdt_annual', 'var') && ~isempty(dhdt_annual)
        disp('Mass conservation enabled - using dhdt_annual for annual mass change conservation.');
    else
        disp('Mass conservation enabled - using h_annual for total ice mass conservation.');
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

    % Determine what to interpolate based on available data
    if exist('dhdt_annual', 'var') && ~isempty(dhdt_annual)
        fprintf('Interpolating dhdt_annual (elevation change) directly...\n');
        data_to_interpolate = dhdt_annual;
        nt = size(dhdt_annual,3);
        data_type = 'elevation change';
    else
        fprintf('Interpolating h_annual (ice thickness) directly...\n');
        data_to_interpolate = h_annual;
        nt = size(dhdt_annual,3)+1;
        data_type = 'ice thickness';
    end

    if size(data_to_interpolate,3) ~= nt
        error('interpolate_altimetry_to_mesh: Length of years must match time dimension of the data_to_interpolate.');
    end

    % Flatten source coordinates
    lat_source_flat = lat_sphere(:);
    long_source_flat = long_sphere(:);

    % Target mesh info
    n_target_vertices = md.mesh.numberofvertices;
    lat_target = md.mesh.lat;
    long_target = md.mesh.long;
    
    % Initialize the spcthickness field with proper dimensions
    md.masstransport.spcthickness = zeros(n_target_vertices + 1, size(dhdt_annual,3)+1);
    fprintf('Initialized spcthickness field: %d vertices x %d time steps\n', n_target_vertices, size(dhdt_annual,3)+1);
    
    % CRITICAL: Define data coverage bounds to prevent global extrapolation
    % Get actual bounds of altimetry data
    data_lat_min = min(lat_sphere(:));
    data_lat_max = max(lat_sphere(:));
    data_lon_min = min(long_sphere(:));
    data_lon_max = max(long_sphere(:));
    
    % Add small buffer (1 degree) to data bounds
    buffer = 1.0;
    data_lat_min = data_lat_min - buffer;
    data_lat_max = data_lat_max + buffer;
    data_lon_min = data_lon_min - buffer;
    data_lon_max = data_lon_max + buffer;
    
    fprintf('Altimetry data coverage: lat %.1f to %.1f, lon %.1f to %.1f\n', ...
            data_lat_min+buffer, data_lat_max-buffer, data_lon_min+buffer, data_lon_max-buffer);
    
    % Create mask for where interpolation should be applied (within data bounds)
    data_coverage_mask = (lat_target >= data_lat_min) & (lat_target <= data_lat_max) & ...
                         (long_target >= data_lon_min) & (long_target <= data_lon_max);
    
    fprintf('Mesh vertices within data coverage: %d (%.1f%% of total)\n', ...
            sum(data_coverage_mask), 100*sum(data_coverage_mask)/length(data_coverage_mask));

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

    % Initialize the interpolated data array
    data_interp_array = zeros(n_target_vertices, nt);

    for t = 1:nt
        % Get current time slice and flatten
        data_slice = data_to_interpolate(:,:,t);
        data_flat = data_slice(:);

        % Identify valid data points for this slice
        valid_mask = ~isnan(lat_source_flat) & ~isnan(long_source_flat) & ~isnan(data_flat);

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
            F = scatteredInterpolant(lat_valid, long_valid, data_valid, 'linear', 'none');
    
            % Interpolate ONLY at mesh points within data coverage
            lat_coverage = lat_target(data_coverage_mask);
            long_coverage = long_target(data_coverage_mask);
            data_coverage = F(lat_coverage, long_coverage);

            % Replace NaN values (from failed interpolation) with zeros
            data_coverage(isnan(data_coverage)) = 0;

            % Additional filtering: only keep values where we actually have nearby data
            % Create a more restrictive mask based on distance to actual data points
            data_mask = ~isnan(data_valid) & abs(data_valid) > 1e-6; % Points with actual data
            if sum(data_mask) > 0
                % Find points that are close to actual data
                lat_data = lat_valid(data_mask);
                long_data = long_valid(data_mask);

                % For each interpolated point, check if it's near data
                max_distance = 0.1; % Maximum distance in degrees (~10 km)
                valid_interp_mask = false(size(data_coverage));

                for i = 1:length(lat_coverage)
                    % Calculate distance to nearest data point
                    distances = sqrt((lat_coverage(i) - lat_data).^2 + (long_coverage(i) - long_data).^2);
                    if min(distances) <= max_distance
                        valid_interp_mask(i) = true;
                    end
                end

                % Zero out interpolated values that are far from data
                data_coverage(~valid_interp_mask) = 0;

                fprintf('    Applied data proximity filter: %d/%d vertices kept (%.1f%%)\n', ...
                        sum(valid_interp_mask), length(data_coverage), 100*sum(valid_interp_mask)/length(data_coverage));
            end

            % Store interpolated values only within coverage area
            data_interp_slice(data_coverage_mask) = data_coverage;
            fprintf('Time step %d: Interpolated at %d vertices (%.1f%% of mesh), %.1f%% have %s > 0\n', ...
                    t, length(data_coverage), 100*length(data_coverage)/n_target_vertices, ...
                    100*sum(abs(data_coverage) > 1e-6)/length(data_coverage), data_type);
        else
            fprintf('  WARNING: No mesh vertices within data coverage area!\n');
        end

        % Update the interpolated data in array
        data_interp_array(:,t) = data_interp_slice;
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

        % Initialize mass conservation report
        mass_conservation_report.original_mass = zeros(nt, 1);
        mass_conservation_report.interpolated_mass = zeros(nt, 1);
        mass_conservation_report.mass_error = zeros(nt, 1);
        mass_conservation_report.correction_factor = zeros(nt, 1);

        for t = 1:nt
            % Original mass calculation based on what was interpolated
            if strcmp(data_type, 'elevation change')
                % We interpolated dhdt_annual - calculate annual mass change
                mass_slice_original = dhdt_annual(:,:,t) * cell_area_original * rhoi * 1e-12; % Convert to Gt
                original_total_mass = sum(mass_slice_original(:), 'omitnan');
                mass_type = 'annual mass change';
            else
                % We interpolated h_annual - calculate total ice mass
                mass_slice_original = h_annual(:,:,t) * cell_area_original * rhoi * 1e-12; % Convert to Gt
                original_total_mass = sum(mass_slice_original(:), 'omitnan');
                mass_type = 'total ice mass';
            end

            % Interpolated mass calculation (from whatever we interpolated)
            data_interp_t =  data_interp_array(:,t);

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

            % Calculate total interpolated mass
            mass_interp_vertices = data_interp_t .* vertex_areas * rhoi * 1e-12; % Convert to Gt
            interpolated_total_mass = sum(mass_interp_vertices, 'omitnan');

            % Debug: Check interpolated data distribution
            non_zero_vertices = sum(abs(data_interp_t) > 1e-6);
            max_data = max(abs(data_interp_t));
            mean_data_nonzero = mean(abs(data_interp_t(abs(data_interp_t) > 1e-6)));
            fprintf('    %s stats: %d vertices with data (%.1f%%), max=%.1f, mean=%.1f\n', ...
                    data_type, non_zero_vertices, 100*non_zero_vertices/n_target_vertices, max_data, mean_data_nonzero);

            % Calculate mass error and correction factor
            mass_error = abs(interpolated_total_mass - original_total_mass);
            if interpolated_total_mass ~= 0
                correction_factor = original_total_mass / interpolated_total_mass;
            else
                correction_factor = 1.0;
            end

            % Debug: Print mass calculation details
            fprintf('  DEBUG - Mass details for time step %d (%s):\n', t, mass_type);
            fprintf('    Original: %.2f Gt\n', original_total_mass);
            fprintf('    Interpolated: %.2f Gt\n', interpolated_total_mass);
            fprintf('    Error: %.2f Gt (%.1f%%)\n', mass_error, 100*mass_error/abs(original_total_mass));
            fprintf('    Correction factor: %.4f\n', correction_factor);

            % Apply mass conservation correction to the interpolated data array
            data_interp_array(:,t) = data_interp_array(:,t) * correction_factor;

            % Calculate corrected mass to verify perfect conservation
            data_corrected_t = data_interp_array(:,t);
            mass_corrected_vertices = data_corrected_t .* vertex_areas * rhoi * 1e-12; % Convert to Gt
            corrected_total_mass = sum(mass_corrected_vertices, 'omitnan');

            % Store statistics
            mass_conservation_report.original_mass(t) = original_total_mass;
            mass_conservation_report.interpolated_mass(t) = interpolated_total_mass;
            mass_conservation_report.corrected_mass(t) = corrected_total_mass;
            mass_conservation_report.mass_error(t) = mass_error;
            mass_conservation_report.correction_factor(t) = correction_factor;
            mass_conservation_report.final_error(t) = abs(corrected_total_mass - original_total_mass);

            fprintf('Time step %d (%s): Original=%.2f Gt, Interpolated=%.2f Gt, Corrected=%.2f Gt\n', ...
                    t, mass_type, original_total_mass, interpolated_total_mass, corrected_total_mass);
            fprintf('              Error before=%.2f Gt (%.1f%%), Error after=%.2e Gt (%.2e%%), Correction=%.4f\n', ...
                    mass_error, 100*mass_error/abs(original_total_mass), ...
                    abs(corrected_total_mass - original_total_mass), ...
                    100*abs(corrected_total_mass - original_total_mass)/abs(original_total_mass), ...
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

        % Determine conservation type for summary
        if exist('dhdt_annual', 'var') && ~isempty(dhdt_annual)
            conservation_type = 'Annual Mass Change Conservation';
        else
            conservation_type = 'Total Ice Mass Conservation';
        end

        fprintf('\n--- %s Summary ---\n', conservation_type);
        fprintf('BEFORE correction - Avg error: %.2f%%, Max error: %.2f%%\n', avg_error_before_percent, max_error_before_percent);
        fprintf('AFTER correction  - Avg error: %.2e%%, Max error: %.2e%%\n', avg_error_after_percent, max_error_after_percent);
        fprintf('Average correction factor: %.4f\n', avg_correction);
        fprintf('=================================\n\n');
    end

    % If we interpolated elevation changes, reconstruct ice thickness time series
    if strcmp(data_type, 'elevation change')
        fprintf('\n=== Reconstructing Ice Thickness Time Series ===\n');
        fprintf('Converting interpolated dhdt_annual to h_annual time series...\n');

        % Create new thickness field by working backwards from assumed final thickness
        h_reconstructed = zeros(n_target_vertices, nt);

        % Set final thickness to a reasonable estimate (using mean from original data if available)
        if ~isempty(h_annual)
            % Use mean final thickness from original data where ice exists
            original_initial = h_annual(:,:,1);
            mean_initial_thickness = mean(original_initial(original_initial > 0), 'omitnan');
            if isnan(mean_initial_thickness) || mean_initial_thickness <= 0
                mean_initial_thickness = 1500; % Default fallback
            end
        else
            mean_initial_thickness = 1500; % Default mean Greenland ice thickness
        end

        fprintf('Using mean initial thickness of %.1f m \n', ...
                mean_initial_thickness);

                 % Work forward through time using dhdt to reconstruct thickness
         h_reconstructed(:, 1) = mean_initial_thickness;
         for t = 1:nt
             h_reconstructed(:, t+1) = h_reconstructed(:, t) + data_interp_array(:,t); % Add dhdt to go forward in time

            % Ensure thickness doesn't go negative where we expect ice
            h_reconstructed(h_reconstructed(:, t) < 0, t) = 0;
        end

        % Update the model structure
        md.masstransport.spcthickness(1:n_target_vertices, :) = h_reconstructed;

        fprintf('Ice thickness time series reconstructed from elevation changes.\n');
        fprintf('Final thickness range: %.1f to %.1f m (mean: %.1f m)\n', ...
                min(h_reconstructed(:, end)), max(h_reconstructed(:, end)), ...
                mean(h_reconstructed(:, end)));

        disp('Interpolation complete - reconstructed h_annual (ice thickness) stored in md.masstransport.spcthickness.');
    else
        disp('Interpolation complete - h_annual (ice thickness) stored in md.masstransport.spcthickness.');
    end

    % Add the time array in the last row of spcthickness
    years_h_annual = [years_dhdt(1)-1, years_dhdt];
    md.masstransport.spcthickness(end, :) = years_h_annual(:)'; % Ensure time is a row vector

    % Return the updated model
    md_updated = md;
    disp('====================================');
end