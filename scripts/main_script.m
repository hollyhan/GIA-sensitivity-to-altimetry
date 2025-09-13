%% Define which steps to run
%steps=[6;7;8;9];
steps = 'refine_mesh_with_gnss';

addpath('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/scripts/mesh');
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
        for i=2:size(h_annual_1,3)
            dhdt = h_annual_1(:,:,i) - h_annual_1(:,:,i-1);
            dhdt1_manual(:,:,i-1) = dhdt;
        end

        for i=1:size(dhdt1_manual,3)
            % compare if dhdt1_manual and dhdt_annual_1 are the same
            if all(dhdt1_manual(:,:,i) == dhdt_annual_1(:,:,i))
                disp('dhdt1_manual and dhdt_annual_1 are the same');
            else
                disp('dhdt1_manual and dhdt_annual_1 are different');
                data = dhdt1_manual(:,:,i) - dhdt_annual_1(:,:,i);
                disp(sum(sum((data))));
                plot_debug = false;
                if plot_debug
                    figure(i)
                    pcolor(lon_vec,lat_vec,dhdt1_manual(:,:,i)-dhdt_annual_1(:,:,i))
                    shading flat;
                    colorbar;
                    title(sprintf('Difference between dhdt1_manual and dhdt_annual_1 for year %d', years_altimetry_1(i)));
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
    % Process glacier mask datasets (returns native high-resolution mask on spherical geographic coordinates)
    disp('=== Processing glacier mask data and apply the mask to the altimetry datasets ===');

    % Find overlapping years between the mask and altimetry data
    yrs_mask = 2022; %1985:2022; % manually defined reliable data years based on the metadata
    % For dhdt, use timestamps starting from the second year (since dhdt = h(t1) - h(t0))
    years_altimetry = [years_altimetry_1(2:end); years_altimetry_2(2:end); years_altimetry_3(2:end); years_altimetry_4(2:end); years_altimetry_5(2:end); years_altimetry_6(2:end); years_altimetry_7(2:end)];
    years_altimetry = unique(years_altimetry);
    yrs_total_overlap = intersect(years_altimetry, yrs_mask);
    data_names = {'measureItsLive', 'DTU', 'Buffalo'};
    data_sets = {dhdt_annual_1, dhdt_annual_2};%, dhdt_annual_3, dhdt_annual_4, dhdt_annual_5, dhdt_annual_6, dhdt_annual_7};

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
`
                    % Mask the dhdt data
                    dhdt_masked{1}(:,:,a) = data_sets{1}(:,:,a).*mask_resampled; % measureItsLive-GEMB
                    dhdt_masked{2}(:,:,a) = data_sets{2}(:,:,a).*mask_resampled; % measureItsLive-GSFC
                    dhdt_masked_years{1}(a) = yrs_total_overlap(n); % Store year for dataset 1
                    dhdt_masked_years{2}(a) = yrs_total_overlap(n); % Store year for dataset 2
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
        save('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results/dhdt_masked.mat', 'dhdt_masked', 'dhdt_masked_years', '-v7.3');
        disp('Saved dhdt_masked.mat with all 7 datasets and their year timetags');
    end

    % Optional plotting (set plot_mask = true to enable)
    plot_mask = false;
    if plot_mask
        for k = 1:length(data_sets)
            if ~isempty(dhdt_masked{k})
                np = 64; % Number of colors
                blue_to_white = [linspace(0,1,np/2)', linspace(0,1,np/2)', ones(np/2,1)];
                white_to_red = [ones(np/2,1), linspace(1,0,np/2)', linspace(1,0,np/2)'];
                custom_cmap = [blue_to_white; white_to_red];

                figure;
                data = dhdt_masked{k}(:,:,1); % Plot first time slice
                pcolor(data);
                shading flat;
                colorbar;
                title(sprintf('Masked dhdt for dataset %d (first time slice)', k));
                xlabel('Longitude', 'FontSize', 14);
                ylabel('Latitude', 'FontSize', 14);
                set(gca, 'FontSize', 14);
                colormap(flip(custom_cmap));
                caxis([-1 1]);
            end
        end

        % Plot the dhdt data and a mask contour 
        figure()
        data = dhdt_annual_1(:,:,end);
        pcolor(long_sphere_1,lat_sphere_1, data)
        shading flat;
        hold on;
        contour(long_sphere_1, lat_sphere_1, mask_resampled, [0.5, 0.5], 'k', 'LineWidth', 2);
        hold off;
        colorbar;
        set(gca, 'FontSize', 14);
        colormap(flip(custom_cmap));
        caxis([-1 1]);
        title('dhdt and a mask contour');
        xlabel('X', 'FontSize', 14);
        ylabel('Y', 'FontSize', 14);

        figure()
        pcolor(long_sphere_1, lat_sphere_1, double(mask_resampled))
        shading flat;
        hold on;
        colorbar;
        set(gca, 'FontSize', 14);
        caxis([-1 1]);
        title('mask field');
        xlabel('X', 'FontSize', 14);
        ylabel('Y', 'FontSize', 14);

        figure()
        A  = dhdt_annual_2(:,:,end);
        M  = mask_resampled;                 % 0/1, same size as A, X_1, Y_1
        D  = A - (A .* M);                   % = -A .* (1 - M)
        % make zero values to NaNs
        D(D == 0) = NaN;
        figure;
        pcolor(X_1, Y_1, D); shading flat; axis equal tight; colorbar;
        hold on; contour(X_1, Y_1, M, [0.5 0.5], 'k', 'LineWidth', 1.5); hold off;
        title('Unmasked minus Masked (EPSG:3413)');
        set(gca, 'FontSize', 14);
        colormap(flip(custom_cmap));
    end

    disp('====================================');
end

if any(steps=='refine_mesh_with_gnss')
    %% GNSS → EPSG:3413, sanity checks, and hierarchical BAMG refinement
    % Assumes you already have:
    %   - lon_gnss, lat_gnss (degrees, WGS84)
    %   - md_regional (ISSM model with an existing 2D mesh in EPSG:3413 meters)
    %     e.g. md_regional = loadmodel(fpath_mesh_model_regional);  % can be anything (e.g., '[]') if bool_mesh_greenland_external is 'false'


    %% ---- 0) Parameters you can tune ----
    r1_km   = 50;         % inner radius (km) at which the mesh is refined to value h1_m
    r2_km   = 100;        % outer radius (km) at which the mesh is refined to value h2_m
    h1_m    = 1e3;        % mesh resolution inside inner radius r1_km, meters
    h2_m    = 5e3;        % mesh resolution between inner and outer radius r1_km and r2_km, meters
    hmax    = 100000;     % mesh resolution beyond radius r2_km, far-field cap (coarse), meters
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
    r2 = r2_km * 1e3;     % 100 km
    h  = hmax * ones(nv,1);             % default far-field size
    h(dmin <= r1) = h1_m;               % ≤ 50 km → 1 km
    mask = dmin > r1 & dmin <= r2;
    h(mask) = h2_m;                     % 50–100 km → 5 km
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

    %% ---- 4) Plot refined mesh with stations ----
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
    pcolor(x_3413_7, y_3413_7, dhdt_annual_7(:,:,end));
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
    pcolor(long_sphere_7, lat_sphere_7, dhdt_annual_7(:,:,end));
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

end
if any(steps==10)
    %addpath('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/scripts/mesh');
    % Load mesh model and initialize models for each altimetry dataset
    %md_regional = loadmodel(fpath_mesh_model_regional);  % can be anything (e.g., '[]') if bool_mesh_greenland_external is 'false'

    % get gnss station coordinates 'lat_gnss' and 'lon_gnss' from running step 1.

    % Assuming the GNSS data sets are already in EPSG:3413 projection, no need to switch the projection.
    % But in case needed, use gdal transform [x,y] = gdaltransform(md.mesh.long,md.mesh.lat,'EPSG:3411','EPSG:3413')

    % Ensure lon in [-180,180), lat in [-90,90]
    lon = mod(lon_gnss+180,360)-180;
    lat = lat_gnss;
    assert(all(isfinite(lon)) && all(isfinite(lat)), 'NaNs/Infs in inputs');
    assert(all(abs(lat) <= 90), 'Lat out of range');

    % Use a PROJ string for WGS84 to force traditional lon,lat order
    [x,y] = gdaltransform(lon, lat, '+proj=longlat +datum=WGS84 +no_defs', 'EPSG:3413');

    fprintf('x range: [%.0f, %.0f] m\n', min(x), max(x));
    fprintf('y range: [%.0f, %.0f] m\n', min(y), max(y));
    % Greenland-typical in EPSG:3413 is roughly:
    % x ~ [-1.5e6, 1.5e6], y ~ [-3.6e6, -2e5]
    [lon2,lat2] = gdaltransform(x, y, 'EPSG:3413', 'EPSG:4326');
    err_lon = max(abs(lon2 - lon));
    err_lat = max(abs(lat2 - lat));
    fprintf('Round-trip max error: lon %.6f°, lat %.6f°\n', err_lon, err_lat);
    figure; scatter(x,y,18,'filled'); axis equal; grid on;
    xlim([-1.5e6 1.5e6]); ylim([-3.6e6 -2e5]);
    xlabel('X (m, EPSG:3413)'); ylabel('Y (m, EPSG:3413)');
    title('GNSS over Greenland — EPSG:3413');

    hmax = 250000; % max edge length in meters
    hmin = 1000; % min edge length in meters
    gradation = 2; % maximum ratio between two adjacent edges
    err = 1;

    md_refined = bamg(md_regional,'hmax',hmax,'hmin',hmin,'gradation',gradation,'splitcorners',1, ...
                      'KeepVertices',0,'RequiredVertices', [x, y]);
    md_refined.mesh.epsg = 3413;
    plotmodel(md_refined,'data','mesh');
    hold on
    plot(x,y,'r.','MarkerSize', 18);
    hold off
    xlabel('X (m, EPSG:3413)'); ylabel('Y (m, EPSG:3413)');
    title('GNSS over Greenland — EPSG:3413');
    legend('Mesh','GNSS');
    set(gca,'FontSize',14);
    axis equal;
end

if any(steps==5)
    addpath('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/scripts/mesh');
    % Load mesh model and initialize models for each altimetry dataset
    md_regional = loadmodel(fpath_mesh_model_regional);  % can be anything (e.g., '[]') if bool_mesh_greenland_external is 'false'

    loading_with_mask = true;
    if loading_with_mask
        % Interpolate with mass conservation enabled and corresponding ice masks
        disp('=== Interpolating with mass conservation enabled and corresponding ice masks ===');
        disp('== Loading measureItsLive-GEMB ==');
        [md1_regional, mass_report_1] = interpolate_altimetry_to_mesh(h_annual_1, lat_sphere_1, long_sphere_1, dhdt_masked_years{1}, md_regional, X_1, Y_1, dhdt_masked{1});

        disp('== Loading measureItsLive-GSFC ==');
        [md2_regional, mass_report_2] = interpolate_altimetry_to_mesh(h_annual_2, lat_sphere_2, long_sphere_2, dhdt_masked_years{2}, md_regional, X_2, Y_2, dhdt_masked{2});

        disp('== Loading DTU2016 ==');
        [md3_regional, mass_report_3] = interpolate_altimetry_to_mesh(h_annual_3, lat_sphere_3, long_sphere_3, dhdt_masked_years{3}, md_regional, X_3, Y_3, dhdt_masked{3});

        disp('== Loading DTU2025 ==');
        [md4_regional, mass_report_4] = interpolate_altimetry_to_mesh(h_annual_4, lat_sphere_4, long_sphere_4, dhdt_masked_years{4}, md_regional, X_4, Y_4, dhdt_masked{4});

        disp('== Loading Buffalo2025-GEMB ==');
        [md5_regional, mass_report_5] = interpolate_altimetry_to_mesh(h_annual_5, lat_sphere_5, long_sphere_5, dhdt_masked_years{5}, md_regional, X_5, Y_5, dhdt_masked{5});

        disp('== Loading Buffalo2025-GSFC ==');
        [md6_regional, mass_report_6] = interpolate_altimetry_to_mesh(h_annual_6, lat_sphere_6, long_sphere_6, dhdt_masked_years{6}, md_regional, X_6, Y_6, dhdt_masked{6});

        disp('== Loading Buffalo2025-IMAU ==');
        [md7_regional, mass_report_7] = interpolate_altimetry_to_mesh(h_annual_7, lat_sphere_7, long_sphere_7, dhdt_masked_years{7}, md_regional, X_7, Y_7, dhdt_masked{7});
    end

    loading_without_mask = false;
    if loading_without_mask
        disp('=== Interpolating without mass conservation and corresponding ice masks ===');
        disp('== Loading measureItsLive-GEMB ==');
        [md1_regional, mass_report_1] = interpolate_altimetry_to_mesh(h_annual_1, lat_sphere_1, long_sphere_1, years_altimetry_1, md_regional, X_1, Y_1, dhdt_annual_1);

        disp('== Loading measureItsLive-GSFC ==');
        [md2_regional, mass_report_2] = interpolate_altimetry_to_mesh(h_annual_2, lat_sphere_2, long_sphere_2, years_altimetry_2, md_regional, X_2, Y_2, dhdt_annual_2);

        disp('== Loading DTU2016 ==');
        [md3_regional, mass_report_3] = interpolate_altimetry_to_mesh(h_annual_3, lat_sphere_3, long_sphere_3, years_altimetry_3, md_regional, X_3, Y_3, dhdt_annual_3);

        disp('== Loading DTU2025 ==');
        [md4_regional, mass_report_4] = interpolate_altimetry_to_mesh(h_annual_4, lat_sphere_4, long_sphere_4, years_altimetry_4, md_regional, X_4, Y_4, dhdt_annual_4);

        disp('== Loading Buffalo2025-GEMB ==');
        [md5_regional, mass_report_5] = interpolate_altimetry_to_mesh(h_annual_5, lat_sphere_5, long_sphere_5, years_altimetry_5, md_regional, X_5, Y_5, dhdt_annual_5);

        disp('== Loading Buffalo2025-GSFC ==');
        [md6_regional, mass_report_6] = interpolate_altimetry_to_mesh(h_annual_6, lat_sphere_6, long_sphere_6, years_altimetry_6, md_regional, X_6, Y_6, dhdt_annual_6);

        disp('== Loading Buffalo2025-IMAU ==');
        [md7_regional, mass_report_7] = interpolate_altimetry_to_mesh(h_annual_7, lat_sphere_7, long_sphere_7, years_altimetry_7, md_regional, X_7, Y_7, dhdt_annual_7);
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

    % plot the mass report (original and interpolated)
    % original mass in solid lines and interpolated mass in dashed lines
    colors = {[0.2157, 0.4941, 0.7216], ...  % Dark Blue for measureItsLive-GEMB
            [0.5294, 0.6667, 0.8627], ...  % Light Blue for measureItsLive-GSFC
            [0.9020, 0.6235, 0.0000], ...  % Dark Orange for DTU2016
            [1.0000, 0.7647, 0.4000], ...  % Light Orange for DTU2025
            [0.4157, 0.2392, 0.6039], ...  % Dark Purple for Buffalo2025-GEMB
            [0.6196, 0.4235, 0.7843], ...  % Medium Purple for Buffalo2025-GSFC
            [0.7725, 0.6392, 0.8706]};     % Light Purple for Buffalo2025-IMAU
    dataset_names = {'measureItsLive-GEMB', 'measureItsLive-GSFC', 'DTU2016', 'DTU2025', 'Buffalo2025-GEMB', 'Buffalo2025-GSFC', 'Buffalo2025-IMAU'};

    figure() % need to debug this for the case of loading with ice mask
    plot(years_altimetry_1(2:end), mass_report_1.original_mass, 'Color', colors{1},'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{1});
    hold on;
    plot(years_altimetry_2(2:end), mass_report_2.original_mass, 'Color', colors{2}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{2});
    plot(years_altimetry_3(2:end), mass_report_3.original_mass, 'Color', colors{3}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{3});
    plot(years_altimetry_4(2:end), mass_report_4.original_mass, 'Color', colors{4}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{4});
    plot(years_altimetry_5(2:end), mass_report_5.original_mass, 'Color', colors{5}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{5});
    plot(years_altimetry_6(2:end), mass_report_6.original_mass, 'Color', colors{6}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{6});
    plot(years_altimetry_7(2:end), mass_report_7.original_mass, 'Color', colors{7}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{7});

    xlabel('Time (year)', 'FontSize', 18);
    ylabel('Mass change (Gt)', 'FontSize', 18);
    title('Annual mass change (Gt/yr) for different altimetry datasets (on native grid)', 'FontSize', 20);
    legend('show', 'Location', 'best', 'FontSize', 16);
    grid on;
    set(gca, 'FontSize', 18);
end

if any(steps==6)
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


if any(steps==7)
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

if any(steps==8)
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

if any(steps==9)
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