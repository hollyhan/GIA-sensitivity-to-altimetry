%% Define which steps to run
steps=[4];

% Load settingsn
run('settings_observation_data.m');

if any(steps==1)
    % Process GNSS data
    [lat_gnss, lon_gnss, data_gnss, err_gnss, time_gnss, R2_gnss] = preprocess_gnss_data(stn_id, fpath_gnss, fname_coord_gnss, n_degree);

    % Array of common time stamps over which to perform model-data comparison
    for n = 1:length(data_gnss)
        if annual_output
            % Find unique years in GNSS data
            common_time{n} = unique(floor(time_gnss{n}));
        else
            % Keep the original gnss data time resolution
            common_time{n} = time_gnss{n};
        end
    end
end

if any(steps==2)
    % Process altimetry datasets without firn correction
    disp('=== Processing altimetry data ===');
    [h_annual_1, dhdt_annual_1, dhdt_monthly_1, years_altimetry_1, lat_sphere_1, long_sphere_1, X_1, Y_1, x_3413_1, y_3413_1] = preprocess_ice_altimetry('measureItsLive-GEMB', false);
    [h_annual_2, dhdt_annual_2, dhdt_monthly_2, years_altimetry_2, lat_sphere_2, long_sphere_2, X_2, Y_2, x_3413_2, y_3413_2] = preprocess_ice_altimetry('measureItsLive-GSFC', false);
    [h_annual_3, dhdt_annual_3, dhdt_monthly_3, years_altimetry_3, lat_sphere_3, long_sphere_3, X_3, Y_3, x_3413_3, y_3413_3] = preprocess_ice_altimetry('DTU2016', false);
    [h_annual_4, dhdt_annual_4, dhdt_monthly_4, years_altimetry_4, lat_sphere_4, long_sphere_4, X_4, Y_4, x_3413_4, y_3413_4] = preprocess_ice_altimetry('DTU2025', false);% DTU data reports-4186.2778 Gt between 2003-2022-12-31 and 4701 Gt if not correcting for firn, Elastic uplift and GIA
    [h_annual_5, dhdt_annual_5, dhdt_monthly_5, years_altimetry_5, lat_sphere_5, long_sphere_5, X_5, Y_5, x_3413_5, y_3413_5] = preprocess_ice_altimetry('Buffalo2025-GEMB', false);
    [h_annual_6, dhdt_annual_6, dhdt_monthly_6, years_altimetry_6, lat_sphere_6, long_sphere_6, X_6, Y_6, x_3413_6, y_3413_6] = preprocess_ice_altimetry('Buffalo2025-GSFC', false);
    [h_annual_7, dhdt_annual_7, dhdt_monthly_7, years_altimetry_7, lat_sphere_7, long_sphere_7, X_7, Y_7, x_3413_7, y_3413_7] = preprocess_ice_altimetry('Buffalo2025-IMAU', false);
    disp('====================================');

    % Debug: compare if dhdt derived from h_annual and dhdt_annual are the same
    debug_dhdt = false;
        if debug_dhdt
        h_to_check = h_annual_7;
        dhdt_to_check = dhdt_annual_7;
        for i=2:size(h_to_check,3)
            dhdt = h_to_check(:,:,i) - h_to_check(:,:,i-1);
            dhdt_manual(:,:,i-1) = dhdt;
        end

        for i=1:size(dhdt_manual,3)
            % compare if dhdt_manual and dhdt_annual are the same
            if all(dhdt_manual(:,:,i) == dhdt_to_check(:,:,i))
                disp('dhdt_manual and dhdt_annual are the same');
            else
                disp('dhdt_manual and dhdt_annual are different');
                data = dhdt_manual(:,:,i) - dhdt_to_check(:,:,i);
                disp(sum(sum((data))));
                plot_debug = false;
                if plot_debug
                    figure(i)
                    pcolor(lon_vec,lat_vec,dhdt_manual(:,:,i)-dhdt_to_check(:,:,i))
                    shading flat;
                    colorbar;
                    title(sprintf('Difference between dhdt_manual and dhdt_annual for year %d', years_altimetry_1(i)));
                    xlabel('Longitude', 'FontSize', 14);
                    ylabel('Latitude', 'FontSize', 14);
                    set(gca, 'FontSize', 14);
                    colormap(flip(custom_cmap))%
                    caxis([-1e-10 1e-10])
                end
            end
        end
    end

end

if any(steps==3)
    plot_mask = false;
    if plot_mask
        np = 64; % Number of colors
        blue_to_white = [linspace(0,1,np/2)', linspace(0,1,np/2)', ones(np/2,1)];
        white_to_red = [ones(np/2,1), linspace(1,0,np/2)', linspace(1,0,np/2)'];
        custom_cmap = [blue_to_white; white_to_red];
    end

    % Process glacier mask datasets (returns native high-resolution mask on spherical geographic coordinates)
    disp('=== Processing glacier mask data and apply the mask to the altimetry datasets ===');

    % Find overlapping years between the mask and altimetry data
    yrs_mask = 1993:2022; % mask data only available from 1993 until 2022
    % For dhdt, use timestamps starting from the second year (since dhdt = h(t1) - h(t0))
    years_altimetry = [years_altimetry_1(2:end); years_altimetry_2(2:end); years_altimetry_3(2:end); years_altimetry_4(2:end); years_altimetry_5(2:end); years_altimetry_6(2:end); years_altimetry_7(2:end)];
    years_altimetry = unique(years_altimetry);
    yrs_total_overlap = intersect(years_altimetry, yrs_mask);
    data_names = {'measureItsLive', 'DTU', 'Buffalo'};
    data_sets = {dhdt_annual_1, dhdt_annual_2, dhdt_annual_3, dhdt_annual_4, dhdt_annual_5, dhdt_annual_6, dhdt_annual_7};

    % Process mask year-by-year to avoid memory issues
    dhdt_masked = cell(length(data_sets), 1);
    dhdt_masked_years = cell(length(data_sets), 1); % Add timetags for each dataset
    
    % Pre-allocate with NaN to track which years have data
    for k = 1:length(data_sets)
        dhdt_masked{k} = NaN(size(data_sets{k},1), size(data_sets{k},2), length(yrs_total_overlap));
        dhdt_masked_years{k} = []; % Initialize empty array for years
    end
    
    % Initialize counters for each dataset group
    a = 1; b = 1; c = 1;
    for n = 2:length(yrs_total_overlap)
        % Extract the native mask data
        if n == 2
            % masks from two consecutive years are needed to derive a union mask
            target_yrs = [yrs_total_overlap(1):yrs_total_overlap(2)];
        else
            target_yrs = yrs_total_overlap(n);
        end

        disp(['== Processing mask for year ', num2str(target_yrs), '==']);
        [mask, mask_year, lat_mask, lon_mask, x_mask, y_mask] = preprocess_glacier_mask(target_yrs);

        % Calculate union mask
        if n == 2
            mask_union = mask(:,:,1) | mask(:,:,2);  % Logical OR for union
            mask_old = mask(:,:,2);
        else
            mask_union = mask | mask_old;  % Logical OR for union
            mask_old = mask;
        end

        for k = 1:length(data_names)
            if strcmp(data_names{k}, 'measureItsLive')
                % check the overlapping years between the mask and the measureItsLive data
                if ismember(yrs_total_overlap(n), years_altimetry_1(2:end))
                    disp('= Resampling glacier mask to the measureItsLive grid =');
                    mask_resampled = resample_mask_to_target_grid_xy(mask_union, x_mask, y_mask, x_3413_2, y_3413_2);
                    mask_resampled = double(mask_resampled); % convert from logical to numeric array

                    % Mask the dhdt data
                    dhdt_masked{1}(:,:,a) = data_sets{1}(:,:,a).*mask_resampled; % measureItsLive-GEMB
                    dhdt_masked{2}(:,:,a) = data_sets{2}(:,:,a).*mask_resampled; % measureItsLive-GSFC
                    dhdt_masked_years{1}(a) = yrs_total_overlap(n); % Store year for dataset 1
                    dhdt_masked_years{2}(a) = yrs_total_overlap(n); % Store year for dataset 2

                    % sanity check
                    if plot_mask
                        figure()
                        pcolor(long_sphere_1, lat_sphere_1, dhdt_masked{1}(:,:,a))
                        shading flat;
                        colorbar;
                        colormap(flip(custom_cmap));
                        caxis([-1 1]);
                        title('dhdt masked{1}');

                        figure() % diff between masked and unmasked
                        pcolor(long_sphere_1, lat_sphere_1, data_sets{1}(:,:,a) - dhdt_masked{1}(:,:,a))
                        shading flat;
                        colormap(flip(custom_cmap));
                        caxis([-1 1]);
                        colorbar;
                        title('unmasked minus masked for measureItsLive-GEMB for a = ', num2str(a));
                    end
                    a = a + 1;
                end
            elseif strcmp(data_names{k}, 'DTU')
                if ismember(yrs_total_overlap(n), years_altimetry_3(2:end))
                    disp('= Resampling glacier mask to the DTU grid =');
                    mask_resampled = resample_mask_to_target_grid_xy(mask_union, x_mask, y_mask, x_3413_3, y_3413_3);
                    mask_resampled = double(mask_resampled); % convert from logical to numeric array
                    dhdt_masked{3}(:,:,b) = data_sets{3}(:,:,b).*mask_resampled; % DTU2016
                    dhdt_masked{4}(:,:,b) = data_sets{4}(:,:,b).*mask_resampled; % DTU2025
                    dhdt_masked_years{3}(b) = yrs_total_overlap(n); % Store year for dataset 3
                    dhdt_masked_years{4}(b) = yrs_total_overlap(n); % Store year for dataset 4

                    % sanity check
                    if plot_mask
                        figure()
                        pcolor(long_sphere_3, lat_sphere_3, dhdt_masked{3}(:,:,b))
                        shading flat;
                        colorbar;
                        colormap(flip(custom_cmap));
                        caxis([-1 1]);
                        title('dhdt masked{3}');

                        figure() % diff between masked and unmasked
                        pcolor(long_sphere_3, lat_sphere_3, data_sets{3}(:,:,b) - dhdt_masked{3}(:,:,b))
                        shading flat;
                        colorbar;
                        colormap(flip(custom_cmap));
                        caxis([-1 1]);
                        title('unmasked minus masked for DTU2025 for b = ', num2str(b));
                    end
                    b = b + 1;
                end
            elseif strcmp(data_names{k}, 'Buffalo')
                if ismember(yrs_total_overlap(n), years_altimetry_5(2:end))
                    disp('= Resampling glacier mask to the Buffalo grid =');
                    mask_resampled = resample_mask_to_target_grid_xy(mask_union, x_mask, y_mask, x_3413_5, y_3413_5);
                    mask_resampled = double(mask_resampled); % convert from logical to numeric array
                    dhdt_masked{5}(:,:,c) = data_sets{5}(:,:,c).*mask_resampled; % Buffalo2025-GEMB
                    dhdt_masked{6}(:,:,c) = data_sets{6}(:,:,c).*mask_resampled; % Buffalo2025-GSFC
                    dhdt_masked{7}(:,:,c) = data_sets{7}(:,:,c).*mask_resampled; % Buffalo2025-IMAU
                    dhdt_masked_years{5}(c) = yrs_total_overlap(n); % Store year for dataset 5
                    dhdt_masked_years{6}(c) = yrs_total_overlap(n); % Store year for dataset 6
                    dhdt_masked_years{7}(c) = yrs_total_overlap(n); % Store year for dataset 7

                    % sanity check
                    if plot_mask
                        figure()
                        pcolor(long_sphere_5, lat_sphere_5, dhdt_masked{5}(:,:,c))
                        shading flat;
                        colorbar;
                        colormap(flip(custom_cmap));
                        caxis([-1 1]);
                        title('dhdt masked{5}');

                        figure() % diff between masked and unmasked
                        pcolor(long_sphere_5, lat_sphere_5, data_sets{5}(:,:,c) - dhdt_masked{5}(:,:,c))
                        shading flat;
                        colorbar;
                        colormap(flip(custom_cmap));
                        caxis([-1 1]);
                        title('unmasked minus masked for Buffalo2025-IMAU for c = ', num2str(c));
                    end
                    c = c + 1;
                end
            end  
        end
    end

    % Remove NaN slices and squeeze the arrays
    disp('=== Cleaning up masked data arrays ===');
    for k = 1:length(data_sets)
        % Find which slices have actual data (not all NaN)
        valid_slices = false(1, size(dhdt_masked{k}, 3));
        for n = 1:size(dhdt_masked{k}, 3)
            if ~all(isnan(dhdt_masked{k}(:,:,n)), 'all')
                valid_slices(n) = true;
            end
        end
     
        % Keep only valid slices
        dhdt_masked{k} = dhdt_masked{k}(:,:,valid_slices);
        dhdt_masked_years{k} = dhdt_masked_years{k}(valid_slices); % Also clean up year array
        disp(['Dataset ', num2str(k), ': kept ', num2str(sum(valid_slices)), ' out of ', num2str(length(valid_slices)), ' time slices']);
    end

    % Attach the first year for h_annual field that will be reconstructed later
    for k = 1:length(dhdt_masked)
        if ~isempty(dhdt_masked{k})
            % Append one earlier timestamp for reconstructing thickness fields
            dhdt_masked_years{k} = [dhdt_masked_years{k}(1) - 1, dhdt_masked_years{k}];
        end
    end

    % DEBUG: Display summary of masked dhdt data
    for k = 1:length(data_sets)
        if ~isempty(dhdt_masked{k})
            nt = size(dhdt_masked{k}, 3);
            disp(['  Years: ', num2str(dhdt_masked_years{k}(1)), ' to ', num2str(dhdt_masked_years{k}(end))]);
            for m = 1:nt
                disp(['  dhdt range (min/max) at time ', num2str(m), ': ', num2str(min(min(dhdt_masked{k}(:,:,m)))), ' to ', num2str(max(max(dhdt_masked{k}(:,:,m)))), ' m/yr']);
            end
        end
    end

    save_data = true;
    if save_data
        save('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results/dhdt_masked_debug.mat', 'dhdt_masked', 'dhdt_masked_years', '-v7.3');
        disp('Saved dhdt_masked with all 7 datasets and their year timetags');
    end

    % Optional plotting (set plot_mask = true to enable)
    if plot_mask
        % Plot the dhdt data and a mask contour 
        figure()
        data = dhdt_annual_7(:,:,end);
        pcolor(long_sphere_7,lat_sphere_7, data)
        shading flat;
        hold on;
        contour(long_sphere_7, lat_sphere_7, mask_resampled, [0.5, 0.5], 'k', 'LineWidth', 2);
        hold off;
        colorbar;
        set(gca, 'FontSize', 14);
        colormap(flip(custom_cmap));
        caxis([-1 1]);
        title('dhdt and a mask contour');
        xlabel('X', 'FontSize', 14);
        ylabel('Y', 'FontSize', 14);

        figure()
        %D = -dhdt_annual_7(:,:,end).*(1 - mask_resampled);
        dhdt_masked{7}(:,:,end) = data_sets{7}(:,:,end).*mask_resampled; % Buffalo2025-IMAU
        D = -dhdt_annual_7(:,:,end) + dhdt_masked{7}(:,:,end);
        % make zero values to NaNs
        D(D == 0) = NaN;
        %contour(long_sphere_7, lat_sphere_7, mask_resampled, [0.5 0.5], '--k', 'LineWidth', 0.2); hold on;
        hold on
        pcolor(long_sphere_7, lat_sphere_7, D);
        shading flat;
        title('Unmasked minus Masked (EPSG:3413)');
        set(gca, 'FontSize', 14);
        colormap(flip(custom_cmap));
        caxis([-1 1]);
        colorbar;
        % plot GNSS stations using lat_gnss and lon_gnss
        plot(lon_gnss, lat_gnss, 'k.', 'MarkerSize', 10);

    end

    disp('====================================');
end

if any(steps==4)
    %% GNSS → EPSG:3413, sanity checks, and hierarchical BAMG refinement
    % Assumes you already have:
    %   - lon_gnss, lat_gnss (degrees, WGS84)
    %   - md_regional (ISSM model with an existing 2D mesh in EPSG:3413 meters)
    %     e.g. md_regional = loadmodel(fpath_mesh_model_regional);  % can be anything (e.g., '[]') if bool_mesh_greenland_external is 'false'

    md_regional = loadmodel(fpath_mesh_model_regional);  % can be anything (e.g., '[]') if bool_mesh_greenland_external is 'false'

    plot_mesh = true;
    if plot_mesh
        np = 64; % Number of colors
        blue_to_white = [linspace(0,1,np/2)', linspace(0,1,np/2)', ones(np/2,1)];
        white_to_red = [ones(np/2,1), linspace(1,0,np/2)', linspace(1,0,np/2)'];
        custom_cmap = [blue_to_white; white_to_red];
    end

    %% ---- 0) Parameters you can tune ----
    r1_km   = 50;         % inner radius (km) at which the mesh is refined to value h1_m
    r2_km   = 150;        % outer radius (km) at which the mesh is refined to value h2_m
    h1_m    = 1e3;        % mesh resolution inside inner radius r1_km, meters
    h2_m    = 5e3;        % mesh resolution between inner and outer radius r1_km and r2_km, meters
    hmax    = 25000;     % mesh resolution beyond radius r2_km, far-field cap (coarse), meters
    hmin    = 1000;       % global minimum edge length, meters. This is a safety clamp to prevent excessively small elements.
    gradation = 2.0;      % bamg: max ratio between neighboring edges. Value 2.0 is a good default ensuring the transition between 1 km and coarser elements is smooth and not too abrupt.

    %% ---- 1) Clean inputs & transform lon/lat → EPSG:3413 (meters) ----
    % Normalize longitudes to [-180, 180), validate lat
    lon = mod(lon_gnss+180,360)-180;
    lat = lat_gnss;
    assert(all(isfinite(lon)) && all(isfinite(lat)), 'NaNs/Infs in inputs');
    assert(all(abs(lat) <= 90), 'Lat out of range');

    % Use PROJ string for 4326 to avoid axis-order ambiguity
    [x, y] = gdaltransform(lon, lat, ...
        '+proj=longlat +datum=WGS84 +no_defs', 'EPSG:3413');

    % Quick numeric sanity check
    fprintf('x range: [%.0f, %.0f] m\n', min(x), max(x));
    fprintf('y range: [%.0f, %.0f] m\n', min(y), max(y));

    % Round-trip check (should be ~meters-level error)
    [lon2, lat2] = gdaltransform(x, y, 'EPSG:3413', 'EPSG:4326');
    err_lon = max(abs(lon2 - lon));
    err_lat = max(abs(lat2 - lat));
    fprintf('Round-trip max error: lon %.6f°, lat %.6f°\n', err_lon, err_lat);

    %% ---- 2) Build per-vertex target size (hVertices) with two refinement rings ----
    % Mesh vertices (meters, EPSG:3413)
    xm = md_regional.mesh.x(:);
    ym = md_regional.mesh.y(:);
    nv = md_regional.mesh.numberofvertices;

    % Station coordinates (EPSG:3413 meters)
    xs = x(:);
    ys = y(:);

    % Compute nearest distance from each mesh vertex to the set of stations
    % Prefer KD-tree if available; otherwise fall back to a robust, memory-safe chunk loop.
    use_kdtree = exist('createns','file') == 2 && exist('knnsearch','file') == 2;
    if use_kdtree
        % Statistics & Machine Learning Toolbox path
        Mdl = createns([xs ys], 'NSMethod', 'kdtree');
        [~, dmin] = knnsearch(Mdl, [xm ym]);   % meters
    else
        % Fallback: chunked min distance to avoid large memory spikes
        fprintf('KD-tree not found; using chunked distance computation...\n');
        chunk = 5000;                      % adjust if you like
        dmin = inf(nv,1);
        for i0 = 1:chunk:nv
            i1 = min(i0+chunk-1, nv);
            XV = xm(i0:i1);
            YV = ym(i0:i1);
            % Compute distances to all stations (vectorized across stations)
            % Dist^2 = (XV - xs')^2 + (YV - ys')^2
            Dx = XV - xs.';
            Dy = YV - ys.';
            D2 = Dx.^2 + Dy.^2;
            dmin(i0:i1) = sqrt(min(D2,[],2));
        end
    end

    % Map distance → target edge length h (meters)
    r1 = r1_km * 1e3;     % 50 km
    r2 = r2_km * 1e3;     % 150 km
    h  = hmax * ones(nv,1);             % default far-field size
    h(dmin <= r1) = h1_m;               % ≤ 50 km → 1 km
    mask = dmin > r1 & dmin <= r2;
    h(mask) = h2_m;                     % 50–150 km → 5 km
    % Global min/max safety clamps
    h = max(h, hmin);
    h = min(h, hmax);

    %% ---- 3) Refine mesh with BAMG, honoring hVertices and locking station vertices ----
    % Ensure unique RequiredVertices
    req = unique([xs ys], 'rows');

    md_refined = bamg(md_regional, ...
        'hVertices', h, ...
        'hmin', hmin, ...
        'hmax', hmax, ...
        'gradation', gradation, ...
        'splitcorners', 1, ...
        'KeepVertices', 1, ...
        'RequiredVertices', req);

    md_refined.mesh.epsg = 3413;

    % Here, also calculate lat and lon values of the refined mesh
    % Define the projection object for EPSG:3413 (WGS 84 / NSIDC Sea Ice Polar Stereographic North)
    proj_info = projcrs(3413);

    % Convert to geographic latitude/longitude (degrees, WGS84)
    [md_refined.mesh.lat, md_refined.mesh.long] = projinv(proj_info, md_refined.mesh.x, md_refined.mesh.y);

    % Lastly, interpolate the geometry and mask fields from the original mesh to the refined mesh
    nVerts_refined = md_refined.mesh.numberofvertices;
    nVerts_coarse = md_regional.mesh.numberofvertices;

    % Define the fields to check
    fields_to_check = {
        'geometry.surface'
        'geometry.base'
        'geometry.thickness'
        'geometry.bed'
        'mask.ocean_levelset'
        'mask.ice_levelset'
    };

    for i = 1:numel(fields_to_check)
        parts = strsplit(fields_to_check{i}, '.');
        value = getfield(md_refined, parts{:});
        field_length = length(value);
        
        if field_length ~= nVerts_refined
            fprintf('md_refined.mesh.numberofvertices is %d, but md_refined.%s has %d vertices. Reinitializing...\n', ...
                    nVerts_refined, fields_to_check{i}, field_length);
            % Reinitialize with zeros
            md_refined = setfield(md_refined, parts{:}, zeros(nVerts_refined, 1));
            
            % --- handle geometry.base separately
            if strcmp(fields_to_check{i}, 'geometry.base') && field_length == nVerts_coarse
                fprintf('→ Interpolating md_regional.%s onto the finer mesh.\n', fields_to_check{i});
                
                Xc = md_regional.mesh.x;
                Yc = md_regional.mesh.y;
                Zc = md_regional.geometry.base;
                Xr = md_refined.mesh.x;
                Yr = md_refined.mesh.y;

                F = scatteredInterpolant(Xc, Yc, Zc, 'natural', 'none');
                md_refined.geometry.base = F(Xr, Yr);
            end
            
            % --- handle mask.ocean_levelset separately
            if strcmp(fields_to_check{i}, 'mask.ocean_levelset') && field_length == nVerts_coarse
                fprintf('→ Interpolating md_regional.%s onto the finer mesh.\n', fields_to_check{i});
                
                Xc = md_regional.mesh.x;
                Yc = md_regional.mesh.y;
                Zc = md_regional.mask.ocean_levelset;
                Xr = md_refined.mesh.x;
                Yr = md_refined.mesh.y;

                F = scatteredInterpolant(Xc, Yc, Zc, 'nearest', 'none');
                md_refined.mask.ocean_levelset = F(Xr, Yr);
            end

            % --- handle mask.ice_levelset separately
            if strcmp(fields_to_check{i}, 'mask.ice_levelset') && field_length == nVerts_coarse
                fprintf('→ Interpolating md_regional.%s onto the finer mesh.\n', fields_to_check{i});
                
                Xc = md_regional.mesh.x;
                Yc = md_regional.mesh.y;
                Zc = md_regional.mask.ice_levelset;
                Xr = md_refined.mesh.x;
                Yr = md_refined.mesh.y;

                F = scatteredInterpolant(Xc, Yc, Zc, 'nearest', 'none');
                md_refined.mask.ice_levelset = F(Xr, Yr);
            end
        end
    end

    %% ---- 4) Plot refined mesh with stations ----
    if plot_mesh
        figure()
        plotmodel(md_refined, 'data', 'mesh');
        hold on; plot(xs, ys, 'b.', 'MarkerSize', 12); hold off;
        axis equal;
        xlabel('X (m, EPSG:3413)'); ylabel('Y (m, EPSG:3413)');
        title('Refined mesh with hierarchical station-centered resolution');
        legend('Mesh','GNSS');
        set(gca,'FontSize',14);

        figure() % plot using x,y coordinates (EPSG:3413)
        % Plot altimetry data first (as background)
        pcolor(x_3413_1, y_3413_1, dhdt_annual_1(:,:,end));
        shading flat;
        hold on; axis equal;
        % Plot mesh lines on top
        triplot(md_refined.mesh.elements, md_refined.mesh.x, md_refined.mesh.y, 'k', 'LineWidth', 0.5);
        % Plot GNSS stations on top
        plot(xs, ys, 'b.', 'MarkerSize', 12); % Use x,y coordinates for GNSS
        xlabel('X (meters, EPSG:3413)'); ylabel('Y (meters, EPSG:3413)');
        title('Refined Mesh + altimetery data field (dhdt) + GNSS stations');
        set(gca,'FontSize',14);
        colorbar;
        colormap(flip(custom_cmap));
        caxis([-1 1]);

        figure() % plot using lat and lon values of the refined mesh
        pcolor(long_sphere_1, lat_sphere_1, dhdt_annual_1(:,:,end));
        shading flat;
        colorbar;
        hold on;
        triplot(md_refined.mesh.elements, md_refined.mesh.long, md_refined.mesh.lat, 'k-');
        plot(lon_gnss, lat_gnss, 'b.', 'MarkerSize', 12);
        colormap(flip(custom_cmap));
        caxis([-1 1]);
        xlabel('Longitude (degrees)'); ylabel('Latitude (degrees)');
        title('Refined Mesh + altimetery data field (dhdt) + GNSS stations');
        set(gca,'FontSize',14);
    end

    %% ---- 5) (Optional) Quick size statistics near stations ----
    fprintf('Refined mesh stats:\n');
    fprintf('  min(hVertices): %.0f m\n', min(h));
    fprintf('  max(hVertices): %.0f m\n', max(h));
    n1 = sum(dmin <= r1);
    n2 = sum(dmin > r1 & dmin <= r2);
    n3 = sum(dmin > r2);
    fprintf('  vertices ≤ %.0f km: %d\n', r1_km, n1);
    fprintf('  vertices in (%.0f–%.0f] km: %d\n', r1_km, r2_km, n2);
    fprintf('  vertices > %.0f km: %d\n', r2_km, n3);

    % Save mesh
    save(fullfile(fpath_mesh_model_regional_refined, 'md_refined.mat'), 'md_regional');
    sprintf('Refined mesh saved in the path: %s',fpath_mesh_model_regional_refined)
end

% Interpolate altimetry data onto the 2D refined regional mesh
if any(steps==5)
    addpath('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/scripts/mesh');

    % For plotting the mass report (original and interpolated)
    % original mass in solid lines and interpolated mass in dashed lines
    colors = {[0.2157, 0.4941, 0.7216], ...  % Dark Blue for measureItsLive-GEMB
            [0.5294, 0.6667, 0.8627], ...  % Light Blue for measureItsLive-GSFC
            [0.9020, 0.6235, 0.0000], ...  % Dark Orange for DTU2016
            [1.0000, 0.7647, 0.4000], ...  % Light Orange for DTU2025
            [0.4157, 0.2392, 0.6039], ...  % Dark Purple for Buffalo2025-GEMB
            [0.6196, 0.4235, 0.7843], ...  % Medium Purple for Buffalo2025-GSFC
            [0.7725, 0.6392, 0.8706]};     % Light Purple for Buffalo2025-IMAU
    dataset_names = {'measureItsLive-GEMB', 'measureItsLive-GSFC', 'DTU2016', 'DTU2025', 'Buffalo2025-GEMB', 'Buffalo2025-GSFC', 'Buffalo2025-IMAU'};

    % Interpolate with mass conservation enabled and corresponding ice masks
    loading_with_mask = true;
    if loading_with_mask
        disp('=== Interpolating with mass conservation enabled and corresponding ice masks ===');
        %disp('== Loading measureItsLive-GEMB ==');
        [md1_regional_with_mask, mass_report_1_dhdt_refined_with_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_1, lat_sphere_1, long_sphere_1, dhdt_masked_years{1}, md_refined, X_1, Y_1, dhdt_masked{1},'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md1_regional_with_mask.mat'), 'md1_regional_with_mask');
        clear md1_regional_with_mask;

        %disp('== Loading measureItsLive-GSFC ==');
        [md2_regional_with_mask, mass_report_2_with_mask] = interpolate_altimetry_to_mesh(h_annual_2, lat_sphere_2, long_sphere_2, dhdt_masked_years{2}, md_regional, X_2, Y_2, dhdt_masked{2});
        save(fullfile(fpath_mesh_model_regional_refined, 'md2_regional_with_mask.mat'), 'md2_regional_with_mask', 'mass_report_2_with_mask');
        clear md2_regional_with_mask;

        disp('== Loading DTU2016 ==');
        [md3_regional_with_mask, mass_report_3_with_mask] = interpolate_altimetry_to_mesh(h_annual_3, lat_sphere_3, long_sphere_3, dhdt_masked_years{3}, md_regional, X_3, Y_3, dhdt_masked{3});
        save(fullfile(fpath_mesh_model_regional_refined, 'md3_regional_with_mask.mat'), 'md3_regional_with_mask', 'mass_report_3_with_mask');
        clear md3_regional_with_mask;

        disp('== Loading DTU2025 ==');
        [md4_regional_with_mask, mass_report_4_with_mask] = interpolate_altimetry_to_mesh(h_annual_4, lat_sphere_4, long_sphere_4, dhdt_masked_years{4}, md_regional, X_4, Y_4, dhdt_masked{4});
        save(fullfile(fpath_mesh_model_regional_refined, 'md4_regional_with_mask.mat'), 'md4_regional_with_mask', 'mass_report_4_with_mask');
        clear md4_regional_with_mask;

        disp('== Loading Buffalo2025-GEMB ==');
        [md5_regional_with_mask, mass_report_5_with_mask] = interpolate_altimetry_to_mesh(h_annual_5, lat_sphere_5, long_sphere_5, dhdt_masked_years{5}, md_regional, X_5, Y_5, dhdt_masked{5});
        save(fullfile(fpath_mesh_model_regional_refined, 'md5_regional_with_mask.mat'), 'md5_regional_with_mask', 'mass_report_5_with_mask');
        clear md5_regional_with_mask;

        disp('== Loading Buffalo2025-GSFC ==');
        [md6_regional_with_mask, mass_report_6_with_mask] = interpolate_altimetry_to_mesh(h_annual_6, lat_sphere_6, long_sphere_6, dhdt_masked_years{6}, md_regional, X_6, Y_6, dhdt_masked{6});
        save(fullfile(fpath_mesh_model_regional_refined, 'md6_regional_with_mask.mat'), 'md6_regional_with_mask', 'mass_report_6_with_mask');
        clear md6_regional_with_mask;

        disp('== Loading Buffalo2025-IMAU ==');
        [md7_regional_with_mask, mass_report_7_with_mask] = interpolate_altimetry_to_mesh(h_annual_7, lat_sphere_7, long_sphere_7, dhdt_masked_years{7}, md_regional, X_7, Y_7, dhdt_masked{7});
        save(fullfile(fpath_mesh_model_regional_refined, 'md7_regional_with_mask.mat'), 'md7_regional_with_mask');
        clear md7_regional_with_mask;

        % save mass report
        save(fullfile(fpath_results_general, 'mass_report_with_mask.mat'), 'mass_report_1_with_mask', 'mass_report_2_with_mask', 'mass_report_3_with_mask', 'mass_report_4_with_mask', 'mass_report_5_with_mask', 'mass_report_6_with_mask', 'mass_report_7_with_mask');

        figure() % need to debug this for the case of loading with ice mask
        plot(years_altimetry_1(2:end), mass_report_1_with_mask.original_mass, 'Color', colors{1},'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{1});
        hold on;
        plot(years_altimetry_2(2:end), mass_report_2_with_mask.original_mass, 'Color', colors{2}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{2});
        plot(years_altimetry_3(2:end), mass_report_3_with_mask.original_mass, 'Color', colors{3}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{3});
        plot(years_altimetry_4(2:end), mass_report_4_with_mask.original_mass, 'Color', colors{4}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{4});
        plot(years_altimetry_5(2:end), mass_report_5_with_mask.original_mass, 'Color', colors{5}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{5});
        plot(years_altimetry_6(2:end), mass_report_6_with_mask.original_mass, 'Color', colors{6}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{6});
        plot(years_altimetry_7(2:end), mass_report_7_with_mask.original_mass, 'Color', colors{7}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{7});

        xlabel('Time (year)', 'FontSize', 18);
        ylabel('Mass change (Gt)', 'FontSize', 18);
        title('Annual mass change (Gt/yr) for different altimetry datasets (on native grid) with ice mask', 'FontSize', 20);
        legend('show', 'Location', 'best', 'FontSize', 16);
        grid on;
        set(gca, 'FontSize', 18);
    end

    loading_without_mask = true;
    if loading_without_mask
        disp('=== Interpolating with mass conservation and without corresponding ice masks ===');
        disp('== Loading measureItsLive-GEMB ==');
        [md1_regional_without_mask, mass_report_1_dhdt_refined_without_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_1, lat_sphere_1, long_sphere_1, years_altimetry_1', md_refined, X_1, Y_1, dhdt_annual_1, 'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md1_regional_without_mask.mat'), 'md1_regional_without_mask');
        clear md1_regional_without_mask;

        disp('== Loading measureItsLive-GSFC ==');
        [md2_regional_without_mask, mass_report_2_without_mask] = interpolate_altimetry_to_mesh(h_annual_2, lat_sphere_2, long_sphere_2, years_altimetry_2', md_regional, X_2, Y_2, dhdt_annual_2);
        save(fullfile(fpath_mesh_model_regional_refined, 'md2_regional_without_mask.mat'), 'md2_regional_without_mask');
        clear md2_regional_without_mask;

        disp('== Loading DTU2016 ==');
        [md3_regional_without_mask, mass_report_3_without_mask] = interpolate_altimetry_to_mesh(h_annual_3, lat_sphere_3, long_sphere_3, years_altimetry_3', md_regional, X_3, Y_3, dhdt_annual_3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md3_regional_without_mask.mat'), 'md3_regional_without_mask');
        clear md3_regional_without_mask;
        
        disp('== Loading DTU2025 ==');
        [md4_regional_without_mask, mass_report_4_without_mask] = interpolate_altimetry_to_mesh(h_annual_4, lat_sphere_4, long_sphere_4, years_altimetry_4', md_regional, X_4, Y_4, dhdt_annual_4);
        save(fullfile(fpath_mesh_model_regional_refined, 'md4_regional_without_mask.mat'), 'md4_regional_without_mask');
        clear md4_regional_without_mask;
        
        disp('== Loading Buffalo2025-GEMB ==');
        [md5_regional_without_mask, mass_report_5_without_mask] = interpolate_altimetry_to_mesh(h_annual_5, lat_sphere_5, long_sphere_5, years_altimetry_5', md_regional, X_5, Y_5, dhdt_annual_5);
        save(fullfile(fpath_mesh_model_regional_refined, 'md5_regional_without_mask.mat'), 'md5_regional_without_mask');
        clear md5_regional_without_mask;
        
        disp('== Loading Buffalo2025-GSFC ==');
        [md6_regional_without_mask, mass_report_6_without_mask] = interpolate_altimetry_to_mesh(h_annual_6, lat_sphere_6, long_sphere_6, years_altimetry_6', md_regional, X_6, Y_6, dhdt_annual_6);
        save(fullfile(fpath_mesh_model_regional_refined, 'md6_regional_without_mask.mat'), 'md6_regional_without_mask');
        clear md6_regional_without_mask;
        
        disp('== Loading Buffalo2025-IMAU ==');
        [md7_regional_without_mask, mass_report_7_without_mask] = interpolate_altimetry_to_mesh(h_annual_7, lat_sphere_7, long_sphere_7, years_altimetry_7', md_regional, X_7, Y_7, dhdt_annual_7);
        save(fullfile(fpath_mesh_model_regional_refined, 'md7_regional_without_mask.mat'), 'md7_regional_without_mask');
        clear md7_regional_without_mask;

        % save mass report
        save(fullfile(fpath_results_general, 'mass_report_without_mask.mat'), 'mass_report_1_without_mask', 'mass_report_2_without_mask', 'mass_report_3_without_mask', 'mass_report_4_without_mask', 'mass_report_5_without_mask', 'mass_report_6_without_mask', 'mass_report_7_without_mask');

        figure() % need to debug this for the case of loading without ice mask
        plot(years_altimetry_1(2:end), mass_report_1_without_mask.original_mass, 'Color', colors{1},'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{1});
        hold on;
        plot(years_altimetry_2(2:end), mass_report_2_without_mask.original_mass, 'Color', colors{2}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{2});
        plot(years_altimetry_3(2:end), mass_report_3_without_mask.original_mass, 'Color', colors{3}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{3});
        plot(years_altimetry_4(2:end), mass_report_4_without_mask.original_mass, 'Color', colors{4}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{4});
        plot(years_altimetry_5(2:end), mass_report_5_without_mask.original_mass, 'Color', colors{5}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{5});
        plot(years_altimetry_6(2:end), mass_report_6_without_mask.original_mass, 'Color', colors{6}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{6});
        plot(years_altimetry_7(2:end), mass_report_7_without_mask.original_mass, 'Color', colors{7}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{7});

        xlabel('Time (year)', 'FontSize', 18);
        ylabel('Mass change (Gt)', 'FontSize', 18);
        title('Annual mass change (Gt/yr) for different altimetry datasets (on native grid) without ice mask', 'FontSize', 20);
        legend('show', 'Location', 'best', 'FontSize', 16);
        grid on;
        set(gca, 'FontSize', 18);
    end
end

if any(steps==6) 
    % create global mesh from regional mesh and initialize model
    loading_with_mask = true;
    if loading_with_mask
        md1_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md1_regional_with_mask.mat'));
        md2_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md2_regional_with_mask.mat'));
        md3_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md3_regional_with_mask.mat'));
        md4_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md4_regional_with_mask.mat'));
        md5_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md5_regional_with_mask.mat'));
        md6_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md6_regional_with_mask.mat'));
        md7_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md7_regional_with_mask.mat'));
    else
        md1_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md1_regional_without_mask.mat'));
        md2_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md2_regional_without_mask.mat'));
        md3_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md3_regional_without_mask.mat'));
        md4_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md4_regional_without_mask.mat'));
        md5_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md5_regional_without_mask.mat'));
        md6_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md6_regional_without_mask.mat'));
        md7_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, 'md7_regional_without_mask.mat'));
    end

    [md1, ~, ~, ~] = createGlobalMesh(md1_regional); % create global mesh from regional mesh
    md1 = initialize_model(md1);

    [md2, ~, ~, ~] = createGlobalMesh(md2_regional); % create global mesh from regional mesh
    md2 = initialize_model(md2);

    [md3, ~, ~, ~] = createGlobalMesh(md3_regional); % create global mesh from regional mesh
    md3 = initialize_model(md3);

    [md4, ~, ~, ~] = createGlobalMesh(md4_regional); % create global mesh from regional mesh
    md4 = initialize_model(md4);

    [md5, ~, ~, ~] = createGlobalMesh(md5_regional); % create global mesh from regional mesh
    md5 = initialize_model(md5);

    [md6, ~, ~, ~] = createGlobalMesh(md6_regional); % create global mesh from regional mesh
    md6 = initialize_model(md6);

    [md7, ~, ~, ~] = createGlobalMesh(md7_regional); % create global mesh from regional mesh
    md7 = initialize_model(md7);
end

if any(steps==7)
    % Calculate GIA using different loading models
    disp('=== Calculating GIA using different loading models ===');

    % Read love numbers from file from path specified in settings_observation_data.m
    load(fpath_love_numbers);

    % Run Green's function method to calculate GIA
    % with ice profile #1
    %[md1, vlm1_VE, vlm1_elastic, hlm1, accm1] = run_gia_greensFunction(md1, ht, lt, kt, false, lat_gnss, lon_gnss);
    md1_solved = run_gia(md1, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md1_solved');

    % with ice profile #2
    %[md2, vlm2_VE, vlm2_elastic, hlm2, accm2] = run_gia_greensFunction(md2, ht, lt, kt, false, lat_gnss, lon_gnss);
    md2_solved = run_gia(md2, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md2_solved');

    % with ice profile #3
    %[md3, vlm3_VE, vlm3_elastic, hlm3, accm3] = run_gia_greensFunction(md3, ht, lt, kt, false, lat_gnss, lon_gnss);
    md3_solved = run_gia(md3, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md3_solved');

    % with ice profile #4
    %[md4, vlm4_VE, vlm4_elastic, hlm4, accm4] = run_gia_greensFunction(md4, ht, lt, kt, false, lat_gnss, lon_gnss);
    md4_solved = run_gia(md4, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md4_solved');

    % with ice profile #5
    %[md5, vlm5_VE, vlm5_elastic, hlm5, accm5] = run_gia_greensFunction(md5, ht, lt, kt, false, lat_gnss, lon_gnss);
    md5_solved = run_gia(md5, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md5_solved');

    % with ice profile #6
    %[md6, vlm6_VE, vlm6_elastic, hlm6, accm6] = run_gia_greensFunction(md6, ht, lt, kt, false, lat_gnss, lon_gnss);
    md6_solved = run_gia(md6, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md6_solved');

    % with ice profile #7
    %[md7, vlm7_VE, vlm7_elastic, hlm7, accm7] = run_gia_greensFunction(md7, ht, lt, kt, false, lat_gnss, lon_gnss);
    md7_solved = run_gia(md7, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md7_solved');

    disp('====================================');
end


if any(steps==8)
    % Compute misfit
    [~, mean1_err_gnss, vlm1_VE_GF_gnss_processed, misfit1, rate1, rate_gnss_fit, y_fit_model1, y_fit_gnss] = ...
    compare_model_to_gnss(lat_gnss, lon_gnss, data_gnss, err_gnss, ...
                        time_gnss, stn_id, md1_solved);

    [~, mean2_err_gnss, vlm2_VE_GF_gnss_processed, misfit2, rate2, rate_gnss_fit, y_fit_model2, y_fit_gnss] = ...
    compare_model_to_gnss(lat_gnss, lon_gnss, data_gnss, err_gnss, ...
                        time_gnss, stn_id, md2_solved);

    [~, mean3_err_gnss, vlm3_VE_GF_gnss_processed, misfit3, rate3, rate_gnss_fit, y_fit_model3, y_fit_gnss] = ...
    compare_model_to_gnss(lat_gnss, lon_gnss, data_gnss, err_gnss, ...
                        time_gnss, stn_id, md3_solved);

    [~, mean4_err_gnss, vlm4_VE_GF_gnss_processed, misfit4, rate4, rate_gnss_fit, y_fit_model4, y_fit_gnss] = ...
    compare_model_to_gnss(lat_gnss, lon_gnss, data_gnss, err_gnss, ...
                        time_gnss, stn_id, md4_solved);

    [~, mean5_err_gnss, vlm5_VE_GF_gnss_processed, misfit5, rate5, rate_gnss_fit, y_fit_model5, y_fit_gnss] = ...
    compare_model_to_gnss(lat_gnss, lon_gnss, data_gnss, err_gnss, ...
                        time_gnss, stn_id, md5_solved);

    [~, mean6_err_gnss, vlm6_VE_GF_gnss_processed, misfit6, rate6, rate_gnss_fit, y_fit_model6, y_fit_gnss] = ...
    compare_model_to_gnss(lat_gnss, lon_gnss, data_gnss, err_gnss, ...
                        time_gnss, stn_id, md6_solved);

    [~, mean7_err_gnss, vlm7_VE_GF_gnss_processed, misfit7, rate7, rate_gnss_fit, y_fit_model7, y_fit_gnss] = ...
    compare_model_to_gnss(lat_gnss, lon_gnss, data_gnss, err_gnss, ...
                        time_gnss, stn_id, md7_solved);
end

if any(steps==9)
    % Plot comparison of all GIA simulations vs GNSS data
    fprintf('\n=== Plotting GIA vs GNSS Comparison ===\n');

    % Define colorblind-friendly colors grouped by altimetry dataset, with different shades for FAC models
    % measureItsLive family (Blue shades)
    colors = {[0.2157, 0.4941, 0.7216], ...  % Dark Blue for measureItsLive-GEMB
              [0.5294, 0.6667, 0.8627], ...  % Light Blue for measureItsLive-GSFC
              [0.9020, 0.6235, 0.0000], ...  % Dark Orange for DTU2016
              [1.0000, 0.7647, 0.4000], ...  % Light Orange for DTU2025
              [0.4157, 0.2392, 0.6039], ...  % Dark Purple for Buffalo2025-GEMB
              [0.6196, 0.4235, 0.7843], ...  % Medium Purple for Buffalo2025-GSFC
              [0.7725, 0.6392, 0.8706]};     % Light Purple for Buffalo2025-IMAU
    dataset_names = {'measureItsLive-GEMB', 'measureItsLive-GSFC', 'DTU2016', 'DTU2025', 'Buffalo2025-GEMB', 'Buffalo2025-GSFC', 'Buffalo2025-IMAU'};

    % Create figure for each GNSS station
    for n = 1:length(data_gnss)
        figure('Position', [100, 100, 1200, 800]);

        % Plot GNSS data points
        plot(time_gnss{n}, data_gnss{n}, 'ko', 'MarkerSize', 8, 'DisplayName', 'GNSS Data', 'LineWidth', 2);
        hold on;

        % Plot GNSS linear fit (should be the same for all comparisons)
        plot(time_gnss{n}, y_fit_gnss{n}, 'k-', 'LineWidth', 3, 'DisplayName', sprintf('GNSS fit (%.2f mm/yr)', rate_gnss_fit(n)));

        % Plot all GIA model fits
        plot(time_gnss{n}, y_fit_model1{n}, 'Color', colors{1}, 'LineWidth', 2, 'DisplayName', sprintf('%s (%.2f mm/yr)', dataset_names{1}, rate1(n)));
        plot(time_gnss{n}, y_fit_model2{n}, 'Color', colors{2}, 'LineWidth', 2, 'DisplayName', sprintf('%s (%.2f mm/yr)', dataset_names{2}, rate2(n)));
        plot(time_gnss{n}, y_fit_model3{n}, 'Color', colors{3}, 'LineWidth', 2, 'DisplayName', sprintf('%s (%.2f mm/yr)', dataset_names{3}, rate3(n)));
        plot(time_gnss{n}, y_fit_model4{n}, 'Color', colors{4}, 'LineWidth', 2, 'DisplayName', sprintf('%s (%.2f mm/yr)', dataset_names{4}, rate4(n)));
        plot(time_gnss{n}, y_fit_model5{n}, 'Color', colors{5}, 'LineWidth', 2, 'DisplayName', sprintf('%s (%.2f mm/yr)', dataset_names{5}, rate5(n)));
        plot(time_gnss{n}, y_fit_model6{n}, 'Color', colors{6}, 'LineWidth', 2, 'DisplayName', sprintf('%s (%.2f mm/yr)', dataset_names{6}, rate6(n)));
        plot(time_gnss{n}, y_fit_model7{n}, 'Color', colors{7}, 'LineWidth', 2, 'DisplayName', sprintf('%s (%.2f mm/yr)', dataset_names{7}, rate7(n)));

        xlabel('Time (year)', 'FontSize', 18);
        ylabel('VLM (mm)', 'FontSize', 18);
        title(sprintf('GIA vs GNSS VLM Rate Comparison at Station %s', stn_id{n}), 'FontSize', 20);
        legend('show', 'Location', 'best', 'FontSize', 16);
        grid on;
        set(gca, 'FontSize', 18);

        % Add text box with statistics
        rates_all = [rate1(n), rate2(n), rate3(n), rate4(n), rate5(n), rate6(n), rate7(n)];
        gnss_rate = rate_gnss_fit(n);
        rate_diff = rates_all - gnss_rate;

        %stats_text = sprintf(['GNSS Rate: %.2f mm/yr\n' ...
        %                    'GIA Rates: %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f mm/yr\n' ...
        %                    'Rate Differences: %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f mm/yr'], ...
        %                    gnss_rate, rates_all, rate_diff);

        %annotation('textbox', [0.02, 0.02, 0.4, 0.15], 'String', stats_text, ...
        %           'EdgeColor', 'black', 'BackgroundColor', 'white', 'FontSize', 10);

        % Save figure with station name in the filename
        saveas(gcf, sprintf('gia_vs_gnss_VLM_rate_at_station_%s.png', stn_id{n}));
    end

    % Create summary figure showing rate comparison across all stations
    figure('Position', [100, 100, 1400, 900]);

    % Prepare data for plotting
    num_stations = length(data_gnss);
    rates_matrix = zeros(num_stations, 7);
    gnss_rates = zeros(num_stations, 1);

    for n = 1:num_stations
        rates_matrix(n, :) = [rate1(n), rate2(n), rate3(n), rate4(n), rate5(n), rate6(n), rate7(n)];
        gnss_rates(n) = rate_gnss_fit(n); % Should be the same for all comparisons
    end

    % Create bar plot
    x_pos = 1:num_stations;
    bar_width = 0.8;

    % Plot GNSS rates as reference dots
    plot(x_pos, gnss_rates, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'black', 'DisplayName', 'GNSS Rate');
    hold on;

    % Plot GIA rates as bars
    for i = 1:7
        bar(x_pos + (i-4)*bar_width/7, rates_matrix(:, i), bar_width/7, ...
            'FaceColor', colors{i}, 'DisplayName', dataset_names{i});
    end

    xlabel('GNSS Station', 'FontSize', 14);
    ylabel('VLM Rate (mm/yr)', 'FontSize', 14);
    title('GIA vs GNSS Rate Comparison Across All Stations', 'FontSize', 16);
    legend('show', 'Location', 'best', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 12);
    xticks(1:num_stations);
    xticklabels(stn_id);
    xtickangle(45);

    % Add statistics summary
    fprintf('\n=== Rate Comparison Summary ===\n');
    fprintf('Station\tGNSS Rate\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', dataset_names{1}, dataset_names{2}, dataset_names{3}, dataset_names{4}, dataset_names{5}, dataset_names{6}, dataset_names{7});
    fprintf('-------\t---------\t---------\t---------\t---------\t---------\t---------\t---------\n');
    for n = 1:num_stations
        fprintf('%s\t%.2f\t\t%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n', ...
                stn_id{n}, gnss_rates(n), rates_matrix(n, :));
    end
end

if any(steps==10)
    % Plot misfits comparison for each station
    figure();
    clf;

    % Prepare data for grouped bar plot - use individual station misfits
    x = 1:length(stn_id);
    %y = [misfit1, misfit2, misfit3, misfit4, misfit5, misfit6, misfit7];  % Each column is a dataset, each row is a station

    % Plot the residual of the model vs the GNSS data (now using vector arrays directly)
    y = [rate_gnss_fit - rate1, rate_gnss_fit - rate2, rate_gnss_fit - rate3, ...
         rate_gnss_fit - rate4, rate_gnss_fit - rate5, rate_gnss_fit - rate6, ...
         rate_gnss_fit - rate7];

    % Create grouped bar plot
    b = bar(x, y, 1);
    hold on;

    % Set colorblind-friendly colors grouped by altimetry dataset, with different shades for FAC models
    % measureItsLive family (Blue shades)
    b(1).FaceColor = [0.2157, 0.4941, 0.7216];  % Dark Blue for measureItsLive-GEMB
    b(2).FaceColor = [0.5294, 0.6667, 0.8627];  % Light Blue for measureItsLive-GSFC
    % DTU family (Orange shades)
    b(3).FaceColor = [0.9020, 0.6235, 0.0000];  % Dark Orange for DTU2016
    b(4).FaceColor = [1.0000, 0.7647, 0.4000];  % Light Orange for DTU2025
    % Buffalo2025 family (Purple shades)
    b(5).FaceColor = [0.4157, 0.2392, 0.6039];  % Dark Purple for Buffalo2025-GEMB
    b(6).FaceColor = [0.6196, 0.4235, 0.7843];  % Medium Purple for Buffalo2025-GSFC  
    b(7).FaceColor = [0.7725, 0.6392, 0.8706];  % Light Purple for Buffalo2025-IMAU

    % Set x-axis labels to station IDs
    set(gca, 'XTick', 1:length(stn_id));
    set(gca, 'XTickLabel', stn_id);

    % Customize the plot to match the reference style
    legend('measureItsLive-GEMB', 'measureItsLive-GSFC', 'DTU2016', 'DTU2025', 'Buffalo2025-GEMB', 'Buffalo2025-GSFC', 'Buffalo2025-IMAU', 'Location', 'northeast');
    xlabel('Station ID');
    xtickangle(45);  % Angled labels for better readability
    ylabel('Residual (mm/yr)');
    title(sprintf('Residual Comparison for different ice loading models (Measured - Model)'));
    set(gca, 'FontSize', 12);
    grid on;
    box on;

    % Set y-axis limits similar to the reference (0 to 50)
    %ylim([0, max(max(y))*1.1]);

    % Save the figure
    saveas(gcf, sprintf('residual_comparison_for_different_ice_loading_models-elastic-with-ice-mask.png'));
end


% Test how much loading differenes will be generated from using masked ice and unmasked ice

% --- Initialize video writer ---
idx_stn = 43;%21; %KUAQ; 43 % KAGA
ylim0 = 69.1;%68.4; %KUAQ %69.1; %KAGA
ylim1 = 69.3;%68.7; %69.3; %KAGA
xlim1 = -49.5;%-33.8; %KUAQ -49.5%KAGA
xlim0 = -50.5; %-32.6; % -50.5;%KAGA
[x_gnss_whole, y_gnss_whole] = ll2xy(lat_gnss, lon_gnss, +1);  % +1 for north polar stereographic (Greenland)
[xmin, ymin] = ll2xy(xlim0, ylim0,+1)
[xmax, ymax] = ll2xy(xlim1, ylim1,+1)

% Differences between masked vs unmasked on the ISSM mesh
nf = 4; 
n0 = 3; % 
t0 = md1_regional_with_mask.masstransport.spcthickness(end,n0);
tf = md1_regional_with_mask.masstransport.spcthickness(end,nf);
diff_h_refined_masked = md1_regional_with_mask.masstransport.spcthickness(1:end-1,nf)-md1_regional_with_mask.masstransport.spcthickness(1:end-1,n0);
diff_h_refined_unmasked =  md1_regional_without_mask.masstransport.spcthickness(1:end-1,nf)-md1_regional_without_mask.masstransport.spcthickness(1:end-1,n0);

plotmodel(md1_regional_with_mask,'data', diff_h_refined_masked,'figure',1,'caxis',([-20 20]), 'xlim',[xmin xmax],'ylim',[ymin ymax])
hold on
title(sprintf('Masked (%d-%d)',t0,tf))
colormap(flip(custom_cmap));
plot(x_gnss, y_gnss, 'k.', 'MarkerSize', 15);

for i = 1:length(stn_id)
    text(x_gnss_whole(i), y_gnss_whole(i), stn_id{i}, ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', 'r', ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center');
end

plotmodel(md1_regional_without_mask,'data', diff_h_refined_unmasked,'figure',2,'caxis',([-20 20]), 'xlim',[xlim0 xlim1],'ylim',[ylim0 ylim1])
hold on
title(sprintf('Unmasked (%d-%d)', md1_regional_without_mask.masstransport.spcthickness(end,n0), md1_regional_without_mask.masstransport.spcthickness(end,nf)))
colormap(flip(custom_cmap));  
plot(x_gnss, y_gnss, 'k.', 'MarkerSize', 15);

for i = 1:length(stn_id)
    text(x_gnss_whole(i), y_gnss_whole(i), stn_id{i}, ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', 'r', ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center');
end

plotmodel(md1_regional_with_mask,'data', diff_h_refined_unmasked-diff_h_refined_masked,'figure',3,'caxis',([-20 20]), 'xlim',[xmin xmax],'ylim',[ymin ymax],'edgecolor','black')
hold on
title(sprintf('Unmasked minus Masked (%d-%d)',t0,tf))
colormap(flip(custom_cmap));  
plot(x_gnss_whole, y_gnss_whole, 'm.', 'MarkerSize', 15);
plot(x_gnss, y_gnss, 'b.', 'MarkerSize', 15);

for i = 1:length(stn_id)
    text(x_gnss_whole(i), y_gnss_whole(i), stn_id{i}, ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', 'r', ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center');
end