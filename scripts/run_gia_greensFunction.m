function [md, vlm_total, vlm_elastic, hlm, accm] = run_gia_greensFunction(md, love_h, love_l, love_k, use_whole_domain, lat_gnss, lon_gnss)
    % Calculate vertical land motion (VLM) using Green's functions
    % Inputs:
    %   md - model structure containing mesh and material properties
    %   love_h - Love number h
    %   love_l - Love number l
    %   love_k - Love number k
    %   use_whole_domain - boolean to use whole domain or only GNSS sites
    %   lat_gnss - latitude of GNSS sites (optional, needed if use_whole_domain is false)
    %   lon_gnss - longitude of GNSS sites (optional, needed if use_whole_domain is false)
    % Outputs:
    %   md - updated model structure
    %   vlm_total - total vertical land motion in meters at each point
    %   vlm_elastic - elastic vertical land motion in meters at each point
    %   hlm - horizontal land motion in meters at each point
    %   accm - acceleration in meters/second^2 at each point

    % Check if GNSS coordinates are provided when needed
    if ~use_whole_domain && nargin < 7
        error('GNSS site coordinates (lat_gnss, lon_gnss) are required when use_whole_domain is false');
    end

    % Load settings
    % -------------
    disp('Reading in the parameterization settings for the run');
    run('settings_gia_parameterization.m');
    run('settings_observation_data.m');

    run_start_year_present = md.masstransport.spcthickness(end,1);
    run_end_year_present = md.masstransport.spcthickness(end,end);
    disp(['Run start year: ', num2str(run_start_year_present)]);
    disp(['Run end year: ', num2str(run_end_year_present)]);
    disp(['Run interval: ', num2str(time_interval_present)]);

    % Get mesh coordinates and areas
    % ------------------------------
    % Get vertex coordinates
    x_vert = md.mesh.x;
    y_vert = md.mesh.y;
    z_vert = md.mesh.z;
    lati_vert = md.mesh.lat;
    long_vert = md.mesh.long;
    
    % Get element coordinates (centroids)
    index_elem = md.mesh.elements;
    x_elem = mean(x_vert(index_elem), 2);
    y_elem = mean(y_vert(index_elem), 2);
    z_elem = mean(z_vert(index_elem), 2);
    
    % Calculate element areas
    a_elem = GetAreas3DTria(index_elem, x_vert, y_vert, z_vert); % m²

    % Calculate distances
    % -------------------
    num_vert = md.mesh.numberofvertices;
    num_elem = md.mesh.numberofelements;
    
    if use_whole_domain == true
        disp('Using the whole Greenland domain for the calculation');
        dx = zeros(num_vert, num_elem);
        % Calculate distance from each vertex to each element centroid
        for jj = 1:num_vert
            x_elem_zero = x_elem - x_vert(jj);
            y_elem_zero = y_elem - y_vert(jj);
            z_elem_zero = z_elem - z_vert(jj);
            dx(jj,:) = sqrt(x_elem_zero.^2 + y_elem_zero.^2 + z_elem_zero.^2);
        end
        n_out = num_vert;
    else   % If GNSS sites are used only
        disp('Using the GNSS sites for the calculation');
        dx = zeros(length(lat_gnss), num_elem);

        % WGS84 ellipsoid parameters
        a = 6378137.0;        % semi-major axis [m]
        f = 1/298.257223563;  % flattening
        b = a*(1-f);          % semi-minor axis [m]
        e2 = 1 - (b/a)^2;     % first eccentricity squared

        % Convert WGS84 lat/lon to spherical lat/lon
        lat_rad = deg2rad(lat_gnss);
        lon_rad = deg2rad(lon_gnss);

        % Calculate geocentric latitude (spherical)
        lat_sph = atan((1-e2) * tan(lat_rad));
        lon_sph = lon_rad;  % Longitude is the same in both systems

        % Convert spherical lat/lon to xyz on sphere
        x_site = r_earth * cos(lat_sph) .* cos(lon_sph);
        y_site = r_earth * cos(lat_sph) .* sin(lon_sph);
        z_site = r_earth * sin(lat_sph);

        % Calculate distance from each element centroid to each GNSS site
        for ii = 1:length(lat_gnss)
            x_site_zero = x_elem - x_site(ii);
            y_site_zero = y_elem - y_site(ii);
            z_site_zero = z_elem - z_site(ii);
            dx(ii,:) = sqrt(x_site_zero.^2 + y_site_zero.^2 + z_site_zero.^2);
        end
        n_out = length(lat_gnss);
    end

    % Calculate scaling factor
    % -------------------------
    rho_e = md.materials.earth_density;  % kg/m³
    rad_e = r_earth; % Earth radius from the 'settings_constants_universal.m')
    scalefactor = 3/(4*pi*rho_e*rad_e^2);  % m/kg (a/m in Farrell 1972 eqn.37)

    % Calculate Green's functions
    % ---------------------------
    % Create distance array for Green's functions (logarithmically spaced)
    dist = exp(linspace(log(100), log(5000*1e3), 1000));  % 100m to 5000km

    % Get the number of time steps in the model
    num_time = length(md.masstransport.spcthickness(end,:)); %length(run_start_year_present:time_interval_present:run_end_year_present);
    nmelt = num_time-1; % Number of loading events
    disp(['Number of time steps in model: ', num2str(num_time)]);

    % Calculate mass changes
    % ----------------------
    % DEBUG: Check ice thickness data before processing
    spcthickness_data = md.masstransport.spcthickness(1:end-1,:);
    fprintf('\n=== DEBUGGING ICE THICKNESS DATA ===\n');
    fprintf('Ice thickness matrix size: %dx%d\n', size(spcthickness_data, 1), size(spcthickness_data, 2));
    fprintf('Total NaN values in ice thickness: %d\n', sum(isnan(spcthickness_data), 'all'));
    fprintf('Total non-zero values: %d\n', sum(abs(spcthickness_data) > 1e-6, 'all'));
    fprintf('Ice thickness range: %.6f to %.6f\n', min(spcthickness_data(:)), max(spcthickness_data(:)));
    fprintf('Sample values from column 1: [%.3f, %.3f, %.3f, ...]\n', spcthickness_data(1,1), spcthickness_data(2,1), spcthickness_data(3,1));
    fprintf('Sample values from column 3: [%.3f, %.3f, %.3f, ...]\n', spcthickness_data(1,3), spcthickness_data(2,3), spcthickness_data(3,3));
    fprintf('=====================================\n\n');

    % Calculate instantaneous mass changes (impulse loads) at each time step
    rho = md.materials.rho_ice;  % kg/m³
    dh_vert = diff(md.masstransport.spcthickness(1:end-1,:), 1, 2);

    % DEBUG: Check the mass changes
    fprintf('=== DEBUGGING MASS CHANGES ===\n');
    fprintf('dh_vert matrix size: %dx%d\n', size(dh_vert, 1), size(dh_vert, 2));
    fprintf('Total NaN values in dh_vert: %d\n', sum(isnan(dh_vert), 'all'));
    fprintf('dh_vert range: %.6f to %.6f\n', min(dh_vert(:)), max(dh_vert(:)));
    fprintf('================================\n\n');
    % Calculate mass changes per element
    % Each column represents instantaneous mass change at that time step
    dm = zeros(num_elem, nmelt);
    for n = 1:nmelt
        dh_vert_slice = dh_vert(:, n);
        dh_elem = mean(dh_vert_slice(index_elem), 2, 'omitnan');  % in meters, ignore NaN
        dm(:, n) = rho * (dh_elem .* a_elem);  % kg (mass change per element)

        % DEBUG: Check for NaN in this time step
        if n <= 3  % Only debug first few time steps to avoid spam
            nan_count_slice = sum(isnan(dh_vert_slice));
            nan_count_elem = sum(isnan(dh_elem));
            nan_count_dm = sum(isnan(dm(:, n)));
            fprintf('Time step %d: dh_vert NaN=%d, dh_elem NaN=%d, dm NaN=%d\n', ...
                    n, nan_count_slice, nan_count_elem, nan_count_dm);
        end
    end

    % DEBUG: Final mass change summary
    total_nan_dm = sum(isnan(dm), 'all');
    total_nonzero_dm = sum(abs(dm) > 1e-6, 'all');
    fprintf('\n=== MASS CHANGE SUMMARY ===\n');
    fprintf('Total NaN in dm: %d\n', total_nan_dm);
    fprintf('Total non-zero in dm: %d\n', total_nonzero_dm);
    fprintf('dm matrix size: %dx%d\n', size(dm, 1), size(dm, 2));
    fprintf('===========================\n\n');

    % Calculate Vertical Land Motion
    % ------------------------------
    % Initialize output arrays
    vlm_elastic = zeros(n_out, nmelt);
    vlm_total = zeros(n_out, nmelt);
    vlm_elastic_accum = zeros(n_out, 1);

    % If elastic only, only need GF with elastic love numbers
    if enable_viscous_deformation_present == false
        % Compute step-response Green's function ONCE for elastic case
        disp('Computing Greens function once for the elastic case');
        [vert, ~, ~] = greens(dist, love_h(:, 1), love_l(:, 1), love_k(:, 1));
        % Interpolate
        gu = interp1(dist, vert, dx);
        gu(isnan(gu)) = 0;

        vlm_accumulated = zeros(n_out, 1);
        % Convolution
        for k = 1:nmelt
            disp(['Calculating elastic only case for loading change step: ', num2str(k)]);
            % Get cumulative load at current time step
            dm_slice = dm(:, k);

            % multiplication
            vlm_elastic(:, k) = scalefactor * (gu * dm_slice);
            vlm_step = scalefactor * (gu * dm_slice);
            vlm_accumulated = vlm_accumulated + vlm_step;
            vlm_elastic(:, k) = vlm_accumulated;
        end
        vlm_total = vlm_elastic;
    else % viscoelastic case
        % Setup Green's function storage
        if use_whole_domain
            % Create directory to save the Green's functions to be memory efficient
            if isunix
                gf_dir = '/tmp/greens_functions';  % Use /tmp on Unix systems for faster I/O
            else
                gf_dir = 'greens_functions';
            end

            if ~exist(gf_dir, 'dir')
                mkdir(gf_dir);
            end
            gu_cell = [];
        else
            % just create a cell array to store the Green's functions
            gu_cell = cell(nmelt, 1);
        end

        % Calculate Green functions first
        for n = 1:nmelt
            disp(['Calculating Green function for loading change step: ', num2str(n)]);
            % Each Green's function represents response to a unit step (Heaviside) load
            [vert, ~, ~] = greens(dist, love_h(:,n), love_l(:,n), love_k(:,n));
            gu = interp1(dist, vert, dx);
            gu(isnan(gu)) = 0;

            if n == 1
                gu_elastic = gu;
            end

            if use_whole_domain
                % Save Green's function to file with compression and v7.3 format
                gf_file = fullfile(gf_dir, sprintf('green_function_%03d.mat', n));
                save(gf_file, 'gu', '-v7.3', '-nocompression');  % v7.3 is faster for large files
            else
                gu_cell{n} = gu;
            end
        end

        % Convolution
        for k = 1:nmelt
            disp(['Calculating viscoelastic response for loading change step: ', num2str(k)]);
            vlm_k = zeros(n_out, 1);  % Reset vlm_k for each time step

            % Calculate total response including viscous effects
            % For each previous load time, calculate its contribution to current response
            for m = 1:k     % For each previous load time
                elapsed_idx = k - m + 1;  % elapsed time index (1 = 0 lag, 2 = 1*dt, ...)

                if use_whole_domain
                    % Load step-response Green's function from file
                    gf_file = fullfile(gf_dir, sprintf('green_function_%03d.mat', elapsed_idx));
                    load(gf_file, 'gu', '-mat');  % Explicitly specify mat format for faster loading
                else
                    gu = gu_cell{elapsed_idx};
                end

                % incremental loading over a timestep
                dm_slice = dm(:, m);

                % simple matrix multiplication
                vlm_timestep = scalefactor * (gu * dm_slice);

                % add the vlm signal from previous loading to complete convolution
                vlm_k = vlm_k + vlm_timestep;
            end

            % --- Accumulate elastic contribution ---
            vlm_elastic_step = scalefactor * (gu_elastic * dm(:, k));
            vlm_elastic_accum = vlm_elastic_accum + vlm_elastic_step;
            vlm_elastic(:, k) = vlm_elastic_accum;

            % --- Total response ---
            vlm_total(:, k) = vlm_k;  % meters
        end
    end

    % Update the model with the calculated vertical land motion
    md.solidearth.lovenumbers.h = love_h;
    md.solidearth.lovenumbers.l = love_l;
    md.solidearth.lovenumbers.k = love_k;

    % Initialize accumulated VLM
    vlm_accumulated = zeros(n_out, 1);

    % Update the model results
    if use_whole_domain
        md.results.TransientSolution = struct('DeltaIceThickness', zeros(num_elem, num_time-1), 'Bed', zeros(num_elem, num_time-1), 'time', zeros(1, num_time-1));
        for n=1:num_time-1
            vlm_accumulated = vlm_accumulated + vlm_total(:,n);
            md.results.TransientSolution(n).DeltaIceThickness = dh_vert(:,n);
            md.results.TransientSolution(n).Bed = md.geometry.bed(:,1) + vlm_accumulated;  % Use initial bed
            md.results.TransientSolution(n).time = run_start_year_present + n*time_interval_present;
        end
    else
        md.results.TransientSolution = struct('DeltaIceThickness', zeros(length(lat_gnss), num_time-1), 'Bed', zeros(length(lat_gnss), num_time-1), 'time', zeros(1, num_time-1));

        % Create interpolant for bed elevation
        F = scatteredInterpolant(md.mesh.x, md.mesh.y, md.geometry.bed(:,1), 'linear');

        for n=1:num_time-1
            vlm_accumulated = vlm_accumulated + vlm_total(:,n);
            md.results.TransientSolution(n).DeltaIceThickness = dh_vert(:,n);
            % Interpolate bed to GNSS sites using scattered data
            bed_gnss = F(x_site, y_site);
            md.results.TransientSolution(n).Bed = bed_gnss + vlm_accumulated;
            md.results.TransientSolution(n).time = run_start_year_present + n*time_interval_present;
        end
    end

    % For now, horizontal and acceleration are not calculated
    hlm = 0;
    accm = 0;
end
