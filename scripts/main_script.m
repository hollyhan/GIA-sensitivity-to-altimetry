% addpaths
addpath('~/Desktop/Code/matlab_functions/arctic-mapping-tools/')
addpath('~/Desktop/Code/matlab_functions/ITS_LIVE/')
addpath('~/Desktop/Data/Altimetry/Greenland/')

%% Define which steps to run
steps=[13];%

% Load settings
run('settings_observation_data.m');
run('settings_gia_parameterization.m');

% Read in external gnss data
if use_berg_et_al
    % update the gnss rates for the 53 stations
    elas_data = read_GNET_Elastic_VLM(fullfile(fpath_gnss_new,'Table_S1_GNET_VLM_Berg_et_al.xlsx'), false);
    %stn_id_berg = cellstr(elas_data.station)';
end

% === Dataset definitions (used across steps) ===
dataset_names = {'MEaSUREs ITS_LIVE-GEMB','MEaSUREs ITS_LIVE-GSFC','DTU2025-RACMO2.3p2', ...
                 'UB-GEMB','UB-GSFC','UB-IMAU'};

% Colors (color-blind-friendly palette)
colors = {[0.2157, 0.4941, 0.7216],[0.5294, 0.6667, 0.8627], ...
          [1.0000, 0.7647, 0.4000],[0.4157, 0.2392, 0.6039],[0.6196, 0.4235, 0.7843], ...
          [0.7725, 0.6392, 0.8706]};

np = 64; % Number of colors
blue_to_white = [linspace(0,1,np/2)', linspace(0,1,np/2)', ones(np/2,1)];
white_to_red = [ones(np/2,1), linspace(1,0,np/2)', linspace(1,0,np/2)'];
custom_cmap = [blue_to_white; white_to_red];

if any(steps==1)
    % Process GNSS data
    [lat_gnss, lon_gnss, data_gnss, err_gnss, time_gnss, R2_gnss] = preprocess_gnss_data(stn_id, fpath_gnss, fname_coord_gnss, n_degree);

    % Array of common time stamps over which to perform model-data comparison
    for n = 1:length(data_gnss)
        if annual_output
            % Find unique years in GNSS data
            if use_berg_et_al
                for n = 1:length(stn_id)
                    name_gnss = string(stn_id{n});
                    match_idx = find(elas_data.station == name_gnss, 1);
                    if ~isempty(match_idx)
                        common_time{n} = (floor(elas_data.tstart(match_idx)): floor(elas_data.tend(match_idx)))';
                    end
                end
            else
                common_time{n} = unique(floor(time_gnss{n}));
            end
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
    [h_annual_3, dhdt_annual_3, dhdt_monthly_3, years_altimetry_3, lat_sphere_3, long_sphere_3, X_3, Y_3, x_3413_3, y_3413_3] = preprocess_ice_altimetry('DTU2025', false);% DTU data reports-4186.2778 Gt between 2003-2022-12-31 and 4701 Gt if not correcting for firn, Elastic uplift and GIA
    [h_annual_4, dhdt_annual_4, dhdt_monthly_4, years_altimetry_4, lat_sphere_4, long_sphere_4, X_4, Y_4, x_3413_4, y_3413_4] = preprocess_ice_altimetry('Buffalo2025-GEMB', false);
    [h_annual_5, dhdt_annual_5, dhdt_monthly_5, years_altimetry_5, lat_sphere_5, long_sphere_5, X_5, Y_5, x_3413_5, y_3413_5] = preprocess_ice_altimetry('Buffalo2025-GSFC', false);
    [h_annual_6, dhdt_annual_6, dhdt_monthly_6, years_altimetry_6, lat_sphere_6, long_sphere_6, X_6, Y_6, x_3413_6, y_3413_6] = preprocess_ice_altimetry('Buffalo2025-IMAU', false);
    disp('====================================');

    % Debug: compare if dhdt derived from h_annual and dhdt_annual are the same
    debug_dhdt = false;
        if debug_dhdt
        h_to_check = h_annual_6;
        dhdt_to_check = dhdt_annual_6;
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
                    title(sprintf('Difference between dhdt_manual and dhdt_annual for year %d', years_altimetry_6(i)));
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

% Bundle per-dataset fields into structured array
load_dhdt_masked_data = true;
if load_dhdt_masked_data
    load('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results/dhdt_masked.mat');
    % If data does not exist, run Step 3
end
datasets = {
    struct('name','MEaSUREs ITS_LIVE-GEMB', 'h',h_annual_1, 'lat',lat_sphere_1, 'lon',long_sphere_1, ...
        'years',years_altimetry_1', 'years_masked', dhdt_masked_years{1}, 'X',X_1, 'Y',Y_1, ...
        'dhdt',dhdt_annual_1, 'dhdt_masked',dhdt_masked{1})
    struct('name','MEaSUREs ITS_LIVE-GSFC', 'h',h_annual_2, 'lat',lat_sphere_2, 'lon',long_sphere_2, ...
        'years',years_altimetry_2', 'years_masked',dhdt_masked_years{2}, 'X',X_2, 'Y',Y_2, ...
        'dhdt',dhdt_annual_2, 'dhdt_masked',dhdt_masked{2})
    struct('name','DTU2025-RACMO2.3p2', 'h',h_annual_3, 'lat',lat_sphere_3, 'lon',long_sphere_3, ...
        'years',years_altimetry_3', 'years_masked',dhdt_masked_years{3}, 'X',X_3, 'Y',Y_3, ...
        'dhdt',dhdt_annual_3, 'dhdt_masked',dhdt_masked{3})
    struct('name','UB-GEMB', 'h',h_annual_4, 'lat',lat_sphere_4, 'lon',long_sphere_4, ...
        'years',years_altimetry_4', 'years_masked',dhdt_masked_years{4}, 'X',X_4, 'Y',Y_4, ...
        'dhdt',dhdt_annual_4, 'dhdt_masked',dhdt_masked{4})
    struct('name','UB-GSFC', 'h',h_annual_5, 'lat',lat_sphere_5, 'lon',long_sphere_5, ...
        'years',years_altimetry_5', 'years_masked',dhdt_masked_years{5}, 'X',X_5, 'Y',Y_5, ...
        'dhdt',dhdt_annual_5, 'dhdt_masked',dhdt_masked{5})
    struct('name','UB-IMAU', 'h',h_annual_6, 'lat',lat_sphere_6, 'lon',long_sphere_6, ...
        'years',years_altimetry_6', 'years_masked',dhdt_masked_years{6}, 'X',X_6, 'Y',Y_6, ...
        'dhdt',dhdt_annual_6, 'dhdt_masked',dhdt_masked{6})
};

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
    years_altimetry = [years_altimetry_1(2:end); years_altimetry_2(2:end); years_altimetry_3(2:end); years_altimetry_4(2:end); years_altimetry_5(2:end); years_altimetry_6(2:end)];
    years_altimetry = unique(years_altimetry);
    yrs_total_overlap = intersect(years_altimetry, yrs_mask);
    data_names = {'measureItsLive', 'DTU', 'Buffalo'};
    data_sets = {dhdt_annual_1, dhdt_annual_2, dhdt_annual_4, dhdt_annual_5, dhdt_annual_6, dhdt_annual_7};

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
                    dhdt_masked{3}(:,:,b) = data_sets{3}(:,:,b).*mask_resampled; % DTU2025
                    dhdt_masked_years{3}(b) = yrs_total_overlap(n); % Store year for dataset 3

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
                    dhdt_masked{4}(:,:,c) = data_sets{4}(:,:,c).*mask_resampled; % Buffalo2025-GEMB
                    dhdt_masked{5}(:,:,c) = data_sets{5}(:,:,c).*mask_resampled; % Buffalo2025-GSFC
                    dhdt_masked{6}(:,:,c) = data_sets{6}(:,:,c).*mask_resampled; % Buffalo2025-IMAU
                    dhdt_masked_years{4}(c) = yrs_total_overlap(n); % Store year for dataset 4
                    dhdt_masked_years{5}(c) = yrs_total_overlap(n); % Store year for dataset 5
                    dhdt_masked_years{6}(c) = yrs_total_overlap(n); % Store year for dataset 6

                    % sanity check
                    if plot_mask
                        figure()
                        pcolor(long_sphere_6, lat_sphere_6, dhdt_masked{6}(:,:,c))
                        shading flat;
                        colorbar;
                        colormap(flip(custom_cmap));
                        caxis([-1 1]);
                        title('dhdt masked{6}');

                        figure() % diff between masked and unmasked
                        pcolor(long_sphere_6, lat_sphere_6, data_sets{6}(:,:,c) - dhdt_masked{6}(:,:,c))
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
        save('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results/dhdt_masked.mat', 'dhdt_masked', 'dhdt_masked_years', '-v7.3');
        disp('Saved dhdt_masked with all 6 datasets and their year timetags');
    end

    % Optional plotting (set plot_mask = true to enable)
    if plot_mask
        % Plot the dhdt data and a mask contour 
        figure()
        data = dhdt_annual_6(:,:,end);
        pcolor(long_sphere_6,lat_sphere_6, data)
        shading flat;
        hold on;
        contour(long_sphere_6, lat_sphere_6, mask_resampled, [0.5, 0.5], 'k', 'LineWidth', 2);
        hold off;
        colorbar;
        set(gca, 'FontSize', 14);
        colormap(flip(custom_cmap));
        caxis([-1 1]);
        title('dhdt and a mask contour');
        xlabel('X', 'FontSize', 14);
        ylabel('Y', 'FontSize', 14);

        figure()
        %D = -dhdt_annual_6(:,:,end).*(1 - mask_resampled);
        dhdt_masked{6}(:,:,end) = data_sets{6}(:,:,end).*mask_resampled; % Buffalo2025-IMAU
        D = -dhdt_annual_6(:,:,end) + dhdt_masked{6}(:,:,end);
        % make zero values to NaNs
        D(D == 0) = NaN;
        %contour(long_sphere_6, lat_sphere_6, mask_resampled, [0.5 0.5], '--k', 'LineWidth', 0.2); hold on;
        hold on
        pcolor(long_sphere_6, lat_sphere_6, D);
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

    plot_mesh = false;
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
        '+proj=longlat +datum=WGS84 +no_defs', 'EPSG:3413'); % if gdaltransform does not work, then simply use the following:[x, y] = projfwd(projcrs(3413), lat, lon);

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
        hold on; plot(xs, ys, 'b.', 'MarkerSize', 16); hold off;
        axis equal;
        xlabel('X (m, EPSG:3413)'); ylabel('Y (m, EPSG:3413)');
        title('Refined mesh with hierarchical station-centered resolution');
        legend('Mesh','GNSS');
        set(gca,'FontSize',14);
        %exportgraphics(gcf, fullfile(fpath_results_figures,'Refined_mesh_with_stations.png'), 'Resolution',300);

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
    save(fullfile(fpath_mesh_model_regional_refined, 'md_refined.mat'), 'md_refined');
    sprintf('Refined mesh saved in the path: %s',fpath_mesh_model_regional_refined)
end

% Interpolate altimetry data onto the 2D refined regional mesh
if any(steps==5)
    addpath('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/scripts/mesh');

    % load the masked data
    load_dhdt_masked_data = false;
    if load_dhdt_masked_data
        load('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results/dhdt_masked.mat');
        % If data does not exist, run Step 3
    end

    % For plotting the mass report (original and interpolated)
    % original mass in solid lines and interpolated mass in dashed lines
    colors = {[0.2157, 0.4941, 0.7216], ...  % Dark Blue for measureItsLive-GEMB
            [0.5294, 0.6667, 0.8627], ...  % Light Blue for measureItsLive-GSFC
            [1.0000, 0.7647, 0.4000], ...  % Light Orange for DTU2025
            [0.4157, 0.2392, 0.6039], ...  % Dark Purple for Buffalo2025-GEMB
            [0.6196, 0.4235, 0.7843], ...  % Medium Purple for Buffalo2025-GSFC
            [0.7725, 0.6392, 0.8706]};     % Light Purple for Buffalo2025-IMAU
    dataset_names = {'measureItsLive-GEMB', 'measureItsLive-GSFC', 'DTU2025', 'Buffalo2025-GEMB', 'Buffalo2025-GSFC', 'Buffalo2025-IMAU'};

    % Interpolate with mass conservation enabled and corresponding ice masks
    loading_with_mask = true;
    if loading_with_mask
        disp('=== Interpolating with mass conservation enabled and corresponding ice masks ===');
        disp('== Loading measureItsLive-GEMB ==');
        [md1_regional_with_mask, mass_report_1_with_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_1, lat_sphere_1, long_sphere_1, dhdt_masked_years{1}, md_refined, X_1, Y_1, dhdt_masked{1},'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md1_regional_with_mask.mat'), 'md1_regional_with_mask');
        clear md1_regional_with_mask;

        disp('== Loading measureItsLive-GSFC ==');
        [md2_regional_with_mask, mass_report_2_with_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_2, lat_sphere_2, long_sphere_2, dhdt_masked_years{2}, md_refined, X_2, Y_2, dhdt_masked{2}, 'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md2_regional_with_mask.mat'), 'md2_regional_with_mask');
        clear md2_regional_with_mask;

        disp('== Loading DTU2025 ==');
        [md3_regional_with_mask, mass_report_3_with_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_3, lat_sphere_3, long_sphere_3, dhdt_masked_years{3}, md_refined, X_4, Y_4, dhdt_masked{3}, 'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md3_regional_with_mask.mat'), 'md3_regional_with_mask');
        clear md3_regional_with_mask;

        disp('== Loading Buffalo2025-GEMB ==');
        [md4_regional_with_mask, mass_report_4_with_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_4, lat_sphere_4, long_sphere_4, dhdt_masked_years{4}, md_refined, X_4, Y_4, dhdt_masked{4}, 'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md4_regional_with_mask.mat'), 'md4_regional_with_mask');
        clear md4_regional_with_mask;

        disp('== Loading Buffalo2025-GSFC ==');
        [md5_regional_with_mask, mass_report_5_with_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_5, lat_sphere_5, long_sphere_5, dhdt_masked_years{5}, md_refined, X_5, Y_5, dhdt_masked{5}, 'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md5_regional_with_mask.mat'), 'md5_regional_with_mask');
        clear md5_regional_with_mask;

        disp('== Loading Buffalo2025-IMAU ==');
        [md6_regional_with_mask, mass_report_6_with_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_6, lat_sphere_6, long_sphere_6, dhdt_masked_years{6}, md_refined, X_6, Y_6, dhdt_masked{6}, 'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md6_regional_with_mask.mat'), 'md6_regional_with_mask');
        clear md6_regional_with_mask;

        % save mass report
        save(fullfile(fpath_results_general, 'mass_report_with_mask.mat'), 'mass_report_1_with_mask', 'mass_report_2_with_mask', 'mass_report_3_with_mask', 'mass_report_4_with_mask', 'mass_report_5_with_mask', 'mass_report_6_with_mask');

        figure(1) % need to debug this for the case of loading with ice mask
        plot(dhdt_masked_years{1}(2:end), mass_report_1_with_mask.mass_src_total, 'Color', colors{1},'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{1});
        hold on;
        plot(dhdt_masked_years{2}(2:end), mass_report_2_with_mask.mass_src_total, 'Color', colors{2}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{2});
        plot(dhdt_masked_years{3}(2:end), mass_report_3_with_mask.mass_src_total, 'Color', colors{3}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{3});
        plot(dhdt_masked_years{4}(2:end), mass_report_4_with_mask.mass_src_total, 'Color', colors{4}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{4});
        plot(dhdt_masked_years{5}(2:end), mass_report_5_with_mask.mass_src_total, 'Color', colors{5}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{5});
        plot(dhdt_masked_years{6}(2:end), mass_report_6_with_mask.mass_src_total, 'Color', colors{6}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{6});

        xlabel('Time (year)', 'FontSize', 18);
        ylabel('Annual mass change (Gt/yr)', 'FontSize', 18);
        title('MB from native altimetry datasets with ice mask', 'FontSize', 20);
        legend('show', 'Location', 'best', 'FontSize', 16);
        grid on;
        set(gca, 'FontSize', 18);
    end

    loading_without_mask = true;
    if loading_without_mask
        disp('=== Interpolating with mass conservation and without corresponding ice masks ===');
        disp('== Loading measureItsLive-GEMB ==');
        [md1_regional_without_mask, mass_report_1_without_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_1, lat_sphere_1, long_sphere_1, years_altimetry_1', md_refined, X_1, Y_1, dhdt_annual_1, 'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md1_regional_without_mask.mat'), 'md1_regional_without_mask');
        clear md1_regional_without_mask;

        disp('== Loading measureItsLive-GSFC ==');
        [md2_regional_without_mask, mass_report_2_without_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_2, lat_sphere_2, long_sphere_2, years_altimetry_2', md_refined, X_2, Y_2, dhdt_annual_2, 'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md2_regional_without_mask.mat'), 'md2_regional_without_mask');
        clear md2_regional_without_mask;

        disp('== Loading DTU2025 ==');
        [md3_regional_without_mask, mass_report_3_without_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_3, lat_sphere_3, long_sphere_3, years_altimetry_3', md_refined, X_3, Y_3, dhdt_annual_3, 'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md3_regional_without_mask.mat'), 'md3_regional_without_mask');
        clear md3_regional_without_mask;
        
        disp('== Loading Buffalo2025-GEMB ==');
        [md4_regional_without_mask, mass_report_4_without_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_4, lat_sphere_4, long_sphere_4, years_altimetry_4', md_refined, X_4, Y_4, dhdt_annual_4, 'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md4_regional_without_mask.mat'), 'md4_regional_without_mask');
        clear md4_regional_without_mask;
        
        disp('== Loading Buffalo2025-GSFC ==');
        [md5_regional_without_mask, mass_report_5_without_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_5, lat_sphere_5, long_sphere_5, years_altimetry_5', md_refined, X_5, Y_5, dhdt_annual_5, 'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md5_regional_without_mask.mat'), 'md5_regional_without_mask');
        clear md5_regional_without_mask;
        
        disp('== Loading Buffalo2025-IMAU ==');
        [md6_regional_without_mask, mass_report_6_without_mask] = interpolate_altimetry_to_mesh_massconservative(h_annual_6, lat_sphere_6, long_sphere_6, years_altimetry_6', md_refined, X_6, Y_6, dhdt_annual_6, 'sigma', 1e3, 'sigma_final', 1e3);
        save(fullfile(fpath_mesh_model_regional_refined, 'md6_regional_without_mask.mat'), 'md6_regional_without_mask');
        clear md6_regional_without_mask;
        
        % save mass report
        save(fullfile(fpath_results_general, 'mass_report_without_mask.mat'), 'mass_report_1_without_mask', 'mass_report_2_without_mask', 'mass_report_3_without_mask', 'mass_report_4_without_mask', 'mass_report_5_without_mask', 'mass_report_6_without_mask');

        figure() % need to debug this for the case of loading without ice mask
        plot(years_altimetry_1(2:end), mass_report_1_without_mask.mass_src_total, 'Color', colors{1},'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{1});
        hold on;
        plot(years_altimetry_2(2:end), mass_report_2_without_mask.mass_src_total, 'Color', colors{2}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{2});
        plot(years_altimetry_3(2:end), mass_report_3_without_mask.mass_src_total, 'Color', colors{3}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{3});
        plot(years_altimetry_4(2:end), mass_report_4_without_mask.mass_src_total, 'Color', colors{4}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{4});
        plot(years_altimetry_5(2:end), mass_report_5_without_mask.mass_src_total, 'Color', colors{5}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{5});
        plot(years_altimetry_6(2:end), mass_report_6_without_mask.mass_src_total, 'Color', colors{6}, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{6});

        xlabel('Time (year)', 'FontSize', 18);
        ylabel('Mass change (Gt)', 'FontSize', 18);
        title('Annual mass change (Gt/yr) for different altimetry datasets (on native grid) without ice mask', 'FontSize', 20);
        legend('show', 'Location', 'best', 'FontSize', 16);
        grid on;
        set(gca, 'FontSize', 18);

        figure()
        plot(dhdt_masked_years{1}(2:end), mass_report_1_with_mask.mass_src_total, 'Color', colors{1},'LineStyle', '-o', 'LineWidth', 2, 'DisplayName', dataset_names{1});
        hold on
        plot(years_altimetry_1(2:end), mass_report_1_without_mask.mass_src_total, 'Color', colors{1},'LineStyle', '-', 'LineWidth', 2, 'DisplayName', dataset_names{1});

    end
end

if any(steps==6) 
    % create global mesh from regional mesh and initialize model
    loading_with_mask = false;
    load_md_regional = true;

    if loading_with_mask
        label = 'with_mask';
    else
        label = 'without_mask';
    end

    if load_md_regional
        md1_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, sprintf('md1_regional_%s.mat', label)));
        md2_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, sprintf('md2_regional_%s.mat', label)));
        md3_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, sprintf('md3_regional_%s.mat', label)));
        md4_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, sprintf('md4_regional_%s.mat', label)));
        md5_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, sprintf('md5_regional_%s.mat', label)));
        md6_regional = loadmodel(fullfile(fpath_mesh_model_regional_refined, sprintf('md6_regional_%s.mat', label)));
    end

    [md1_global, ~, ~, ~] = createGlobalMesh(md1_regional); % create global mesh from regional mesh
    md1_global = initialize_model(md1_global);
    save(fullfile(fpath_mesh_model_regional_refined, sprintf('md1_global_%s.mat', label)), 'md1_global');
    clear md1_regional

    [md2_global, ~, ~, ~] = createGlobalMesh(md2_regional); % create global mesh from regional mesh
    md2_global = initialize_model(md2_global);
    save(fullfile(fpath_mesh_model_regional_refined, sprintf('md2_global_%s.mat', label)), 'md2_global');
    clear md2_regional

    [md3_global, ~, ~, ~] = createGlobalMesh(md3_regional); % create global mesh from regional mesh
    md3_global = initialize_model(md3_global);
    save(fullfile(fpath_mesh_model_regional_refined, sprintf('md3_global_%s.mat', label)), 'md3_global');
    clear md3_regional

    [md4_global, ~, ~, ~] = createGlobalMesh(md4_regional); % create global mesh from regional mesh
    md4_global = initialize_model(md4_global);
    save(fullfile(fpath_mesh_model_regional_refined, sprintf('md4_global_%s.mat', label)), 'md4_global');
    clear md4_regional

    [md5_global, ~, ~, ~] = createGlobalMesh(md5_regional); % create global mesh from regional mesh
    md5_global = initialize_model(md5_global);
    save(fullfile(fpath_mesh_model_regional_refined, sprintf('md5_global_%s.mat', label)), 'md5_global');
    clear md5_regional

    [md6_global, ~, ~, ~] = createGlobalMesh(md6_regional); % create global mesh from regional mesh
    md6_global = initialize_model(md6_global);
    save(fullfile(fpath_mesh_model_regional_refined, sprintf('md6_global_%s.mat', label)), 'md6_global');
    clear md6_regional
end

if any(steps==7)
    % Calculate GIA using different loading models
    disp('=== Calculating GIA using different loading models ===');

    loading_with_mask = true;
    load_md_global = true;
    save_md = true;

    if loading_with_mask
        label = 'with_mask';
    else
        label = 'without_mask';
    end

    if load_md_global
        mds = cell(numel(datasets), 1);
        for k = 1:numel(datasets)
            ds = datasets{k};
            fname = fullfile(fpath_mesh_model_regional_refined, sprintf('md%d_global_%s.mat', k, label));
            disp(sprintf('loading md data for %s', ds.name));
            mds{k} = loadmodel(fname);
        end
    end

    % Read love numbers from file from path specified in settings_observation_data.m
    load(fpath_love_numbers);

    % Run Green's function method to calculate GIA
    mds_solved = cell(numel(datasets), 1);
    vlm_total = cell(numel(datasets), 1);
    vlm_elastic = cell(numel(datasets), 1);
    hlm = cell(numel(datasets), 1);
    accm = cell(numel(datasets), 1);

    for k = 1:numel(datasets)
        [mds_solved{k}, vlm_total{k}, vlm_elastic{k}, hlm{k}, accm{k}] = run_gia_greensFunction(mds{k}, ht, lt, kt, false, lat_gnss, lon_gnss);
        %new_name = sprintf('md%d_solved', k)
        %mds_solved{k} = run_gia(mds{k}, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, new_name);

        if save_md
            md_save = mds_solved{k};
            vlm_save = vlm_total{k};
            disp(sprintf('Saving solved model object to file: %s', fullfile(fpath_results_general, sprintf('md%d_solved_%s_%s.mat', k, rheology_choice, label))));
            save(fullfile(fpath_results_general, sprintf('md%d_solved_%s_%s.mat', k, rheology_choice, label)), 'md_save');
            save(fullfile(fpath_results_general, sprintf('vlm%d_total_%s_%s.mat', k, rheology_choice, label)), 'vlm_save');
        end
    end
    disp('====================================');
end

if any(steps==8)

    load_md_solved = true;
    if load_md_solved
        mds_solved = cell(numel(datasets), 1);
        for k = 1:numel(datasets)
            fname = fullfile(fpath_results_general, sprintf('md%d_solved_%s_%s.mat', k, rheology_choice, label));
            disp(sprintf('loading md_solved data %s for %s', fname, datasets{k}.name));
            mds_solved{k} = loadmodel(fname);
        end
        vlm_total = cell(numel(datasets), 1);
        for k = 1:numel(datasets)
            fname = fullfile(fpath_results_general, sprintf('vlm%d_total_%s_%s.mat', k, rheology_choice, label));
            disp(sprintf('loading vlm_total data %s for %s', fname, datasets{k}.name));
            vlm_total{k} = load(fname).vlm_save;
        end
    end

    % Compute misfits between all model results and GNSS data
    vlm_VE_GF_gnss  = cell(numel(datasets), 1);
    misfit          = cell(numel(datasets), 1);
    rate_model      = cell(numel(datasets), 1);
    y_fit_model     = cell(numel(datasets), 1);
    y_fit_gnss      = cell(numel(datasets), 1);

    for k = 1:numel(datasets)
        disp('------------------------------------')
        disp(sprintf('Comparing dataset %s to GNSS data...\n', datasets{k}.name));
        if use_berg_et_al
            [~, mean_err_gnss, vlm_VE_GF_gnss{k}, misfit{k}, ...
                rate_model{k}, rate_gnss_fit, y_fit_model{k}, ~] = ...
                compare_model_to_gnss(lat_gnss, lon_gnss, data_gnss, err_gnss, ... % err_gnss might need to be revised to reflect Danjals err
                                    common_time, stn_id, mds_solved{k}, vlm_total{k});
        else
            [~, mean_err_gnss, vlm_VE_GF_gnss{k}, misfit{k}, ...
                rate_model{k}, rate_gnss_fit, y_fit_model{k}, y_fit_gnss{k}] = ...
                compare_model_to_gnss(lat_gnss, lon_gnss, data_gnss, err_gnss, ...
                                    common_time, stn_id, mds_solved{k}, vlm_total{k});
        end
    end
end

if any(steps==9)
    save_fig = false;
    % Plot comparison of all GIA simulations vs GNSS data (Fig. 4 and Fig. Supp 2)
    fprintf('\n=== Plotting GIA vs GNSS Comparison ===\n');
    % Define colorblind-friendly colors grouped by altimetry dataset, with different shades for FAC models
    % measureItsLive family (Blue shades)
    colors = {[0.2157, 0.4941, 0.7216], ...  % Dark Blue for measureItsLive-GEMB
              [0.5294, 0.6667, 0.8627], ...  % Light Blue for measureItsLive-GSFC
              [1.0000, 0.7647, 0.4000], ...  % Light Orange for DTU2025
              [0.4157, 0.2392, 0.6039], ...  % Dark Purple for Buffalo2025-GEMB
              [0.6196, 0.4235, 0.7843], ...  % Medium Purple for Buffalo2025-GSFC
              [0.7725, 0.6392, 0.8706]};     % Light Purple for Buffalo2025-IMAU

    num_stations = length(data_gnss);
    num_datasets = numel(datasets);

    % Create figure for each GNSS station (this can be done only using the gnss data processed in Step 1, not using Berg et al.s datasets)
    if ~use_berg_et_al
        for n = 1:num_stations
            figure('Position', [100, 100, 1200, 800]);

            % Plot GNSS data points
            plot(time_gnss{n}, data_gnss{n}, 'ko', 'MarkerSize', 8, 'DisplayName', 'GNSS Data', 'LineWidth', 2);
            hold on;

            % Plot GNSS linear fit (y_fit_gnss should be the same, so just pick one)
            plot(time_gnss{n}, y_fit_gnss{1}{n}, 'k-', 'LineWidth', 3, 'DisplayName', sprintf('GNSS fit (%.2f mm/yr)', rate_gnss_fit(n)));

            % Plot all GIA model fits
            for k = 1:num_datasets
                plot(time_gnss{n}, y_fit_model{k}{n}, 'Color', colors{k}, 'LineWidth', 2, ...
                    'DisplayName', sprintf('%s (%.2f mm/yr)', ...
                    datasets{k}.name, rate_model{k}(n)));
            end

            xlabel('Time (year)', 'FontSize', 18);
            ylabel('VLM (mm)', 'FontSize', 18);
            title(sprintf('GIA vs GNSS VLM Rate Comparison at Station %s (%s)', stn_id{n}, label), 'FontSize', 20);
            legend('show', 'Location', 'best', 'FontSize', 16);
            grid on;
            set(gca, 'FontSize', 18);

            % Save figure with station name in the filename
            if save_fig
                saveas(gcf, fullfile(fpath_results_figures, sprintf('gia_vs_gnss_VLM_rate_at_station_%s_%s_%s.png', stn_id{n}, rheology_choice, label)));
            end
        end
    end

    % Create summary figure showing rate comparison across all stations
    % Prepare data for plotting
    if use_berg_et_al
        gnss_rates_berg = NaN(num_stations,1);
        gnss_rates_err_berg = NaN(num_stations,1);
        for i = 1:num_stations
            name_gnss = string(stn_id{i});
            match_idx = find(elas_data.station == name_gnss, 1);
            if ~isempty(match_idx)
                gnss_rates_berg(i) = elas_data.Uobserved(match_idx);
                gnss_rates_err_berg(i) = elas_data.Uobserved_sigma(match_idx);
            end
        end
        rate_gnss = gnss_rates_berg;
    else
        rate_gnss = rate_gnss_fit;
    end

    rates_matrix = zeros(num_stations, num_datasets);
    for k = 1:num_datasets
        rates_matrix(:, k) = rate_model{k};
    end

    % Create bar plot (Fig. Supp 2)
    x_pos = 1:num_stations;
    bar_width = 0.8;

    figure('Color','w','Position',[100 100 1700 600]);
    % Plot GIA rates as bars
    for k = 1:num_datasets
        bar(x_pos + (k - (num_datasets+1)/2) * bar_width/num_datasets, ...
            rates_matrix(:, k), bar_width/num_datasets, ...
            'FaceColor', colors{k}, 'DisplayName', datasets{k}.name);
        hold on
    end

    % Plot GNSS rates as reference dots
    % Figure S5
    errorbar(x_pos, rate_gnss, gnss_rates_err_berg, ...
         'o', 'Color','k','MarkerFaceColor','none','MarkerSize',5.5, 'LineWidth',1.2,'DisplayName', 'Observed GNSS rate (Berg et al., 2024)');
    %plot(x_pos, rate_gnss, 'ko', 'MarkerSize', 5.5, 'MarkerFaceColor', 'black', 'DisplayName', 'Observed GNSS rate');


    % --- Labels & layout ---
    xlabel('Station ID', 'FontSize', 14);
    ylabel('VLM Rate (mm/yr)', 'FontSize', 14);
    title('Modeld elastic and GNSS-observed VLM rates at Greenland GNSS stations', 'FontSize', 16);
    legend('show', 'Location', 'northwest', 'FontSize', 14, 'Interpreter','none');
    grid on;
    set(gca, 'FontSize', 14);
    xticks(x_pos);
    xticklabels(stn_id);
    xtickangle(45);

    % Save the figure
    if save_fig
        saveas(gcf, fullfile(fpath_results_figures, sprintf('FigureS5-Elastic vs GNSS Rate Comparison Across All Stations_%s_%s.png', rheology_choice, label)));
    end

    % Add statistics summary
    fprintf('\n=== Rate Comparison Summary ===\n');
    fprintf('Station\tGNSS\t');
    dataset_names = cellfun(@(d) d.name, datasets, 'UniformOutput', false);
    fprintf('%s\t', dataset_names{:});
    fprintf('\n--------------------------------------------------------------\n');
    for n = 1:num_stations
        fprintf('%s\t%.2f\t', stn_id{n}, gnss_rates_berg(n));
        fprintf('%.2f\t', rates_matrix(n, :));
        fprintf('\n');
    end
end

if any(steps == 10)
    fprintf('\n=== Plotting GIA–GNSS Residual Comparison ===\n');
    % Figure S4.
    save_fig = true;

    num_datasets = numel(datasets);
    num_stations = numel(stn_id);

    %% --- Compute residuals (GNSS − Model) ---
    if use_berg_et_al
        gnss_rates_berg = NaN(1,num_stations);
        gnss_rates_err_berg = NaN(1,num_stations);
        for i = 1:num_stations
            name_gnss = string(stn_id{i});
            match_idx = find(elas_data.station == name_gnss, 1);
            if ~isempty(match_idx)
                gnss_rates_berg(i) = elas_data.Uobserved(match_idx);
                gnss_rates_err_berg(i) = elas_data.Uobserved_sigma(match_idx);
            end
        end
        rate_gnss = gnss_rates_berg;
    else
        rate_gnss = rate_gnss_fit;
    end

    y = zeros(num_stations, num_datasets);
    for k = 1:num_datasets
        y(:, k) = rate_gnss - rate_model{k};
    end

    %% --- Read & apply the paleo and contemporary loading corrections ---
    % Data from Parviz and Glenn for the paleo signals
    gia_data = read_VLM_sites(fpath_gia, false);
    gia1D_total_mean = gia_data.VLM1D_mean; % this includes total signal (DG+LIA+PG)
    gia3D_total_mean = gia_data.VLM3D_mean;
    gia1D_total_mean_sigma = gia_data.VLM1D_sigma;
    gia3D_total_mean_sigma = gia_data.VLM3D_sigma;
    gia_DG3D_mean = gia_data.DG; % mean of 3D Deglacial-GIA
    gia_DG3D_sigma = gia_data.DG_sigma;
    gia_PG = gia_data.PG;
    gia_LIA_mean = gia_data.LIA_mean;
    gia_LIA_sigma = gia_data.LIA_sigma;

    y_corr_1D_total = NaN(size(y));
    y_corr_3D_total = NaN(size(y));
    y_corr_3D_DG_PG = NaN(size(y));
    applied_gia_1D_total = NaN(num_stations,1);
    applied_gia_3D_total = NaN(num_stations,1);
    applied_gia_3D_DG_PG = NaN(num_stations,1);
    gia_val_LIA = NaN(num_stations,1);

    for i = 1:num_stations
        name_gnss = string(stn_id{i});
        match_idx = find(gia_data.name == name_gnss, 1);
        if ~isempty(match_idx)
            gia_val_1D_total = gia1D_total_mean(match_idx);
            gia_val_3D_total = gia3D_total_mean(match_idx);
            gia_val_3D_DG_PG = gia_DG3D_mean(match_idx) + gia_PG(match_idx);
            gia_val_LIA(i) = gia_data.LIA_from_residual(match_idx); % Note: this can change depending on Parvizs ansnwer to my question
            y_corr_1D_total(i,:) = y(i,:) - gia_val_1D_total;
            y_corr_3D_total(i,:) = y(i,:) - gia_val_3D_total;
            y_corr_3D_DG_PG(i,:) = y(i,:) - gia_val_3D_DG_PG;
            applied_gia_1D_total(i) = gia_val_1D_total;
            applied_gia_3D_total(i) = gia_val_3D_total;
            applied_gia_3D_DG_PG(i) = gia_val_3D_DG_PG;
        else
            warning('⚠️ No GIA match found for station %s', name_gnss);
        end
    end

    nMatched = sum(~isnan(applied_gia_1D_total));
    fprintf('\n✅ GIA corrections applied for %d of %d stations.\n', nMatched, num_stations);

    % Use the data from Berg et al. (2024) for the contemporary Perpheral Glacier Loading
    % Apply the elastic corrections to the residuals corrected for the total 3D GIA signal
    y_corr_total = NaN(size(y));
    applied_elas_Can = NaN(num_stations,1);
    comparison_elas = NaN(size(y));
    for i = 1:num_stations
        name_gnss = string(stn_id{i});
        match_idx = find(elas_data.station == name_gnss, 1);
        if ~isempty(match_idx)
            fprintf('Applying corrections to station %s \n', name_gnss)
            elas_Can = elas_data.Uelastic_CanPG(match_idx); % only consider the Canadian component because the Greenland component is already included
            elas_Can_sigma = elas_data.Uelastic_CanPG_sigma(match_idx);
            applied_elas_Can(i) = elas_Can;
            y_corr_total(i,:) = y(i,:) - applied_gia_3D_total(i) - applied_elas_Can(i);
            comparison_elas(i,:) = rates_matrix(i,:) - elas_data.Uelastic_GrIS(match_idx);
            gnss_rates_berg(i) = elas_data.Uobserved(match_idx);
            gnss_rates_err_berg(i) = elas_data.Uobserved_sigma(match_idx);
        else
            warning('⚠️ No Elastic GrCan match found for station %s', name_gnss);
        end
    end

    nMatched = sum(~isnan(applied_elas_Can));
    fprintf('\n✅ Elastic Greenland and Canadian PG corrections applied for %d of %d stations.\n', nMatched, num_stations);

    %% --- Plot raw residuals (Figure S4a)---
    figure('Color', 'w', 'Position', [100 100 1700 600]);
    x = 1:num_stations;
    b = bar(x, y, 1);
    hold on;

    colors = { [0.2157,0.4941,0.7216], [0.5294,0.6667,0.8627], ...
               [1.0000,0.7647,0.4000], ...
               [0.4157,0.2392,0.6039], [0.6196,0.4235,0.7843], ...
               [0.7725,0.6392,0.8706] };
    for k = 1:num_datasets
        b(k).FaceColor = colors{k};
    end

    set(gca,'XTick',x,'XTickLabel',stn_id,'FontSize',14);
    xtickangle(45);
    xlabel('Station ID'); ylabel('Residual (mm/yr)');
    title('Residuals (GNSS minus Model)');
    legend(cellfun(@(d)d.name,datasets,'UniformOutput',false),'Location','northwest','FontSize', 14,'Interpreter','none');
    grid on; box on;
    ylim([-7 13]);
    if save_fig
        saveas(gcf, fullfile(fpath_results_figures, sprintf('Residuals_Raw_%s.png', label)));
    end

    %% --- Quick check summary for a single altimetry product (choose the first dataset)---
    nDigits = 2;
    T = table(stn_id(:), ...
            round(y(:,1), nDigits), ... % Raw residual
            round(applied_gia_1D_total, nDigits), ...
            round(applied_gia_3D_DG_PG, nDigits), ...
            round(applied_gia_3D_total, nDigits), ...
            round(applied_elas_Can, nDigits), ...
            round(y_corr_1D_total(:,1), nDigits), ...
            round(y_corr_3D_DG_PG(:,1), nDigits), ...
            round(y_corr_3D_total(:,1), nDigits), ...
            round(y_corr_total(:,1), nDigits), ...
            'VariableNames', {'Station','Raw residual', 'GIA_1D_total','GIA_3D_DG_PG','GIA_3D_total', ...
                                'contemporary PG (Can)','Corrected_1D_total', ...
                                'Corrected_3D_DG_PG','Corrected_3D_total','Corrected_total'});
    disp(sprintf('Displaying applied signals and corrected residuals (mm/yr for %s)', datasets{1}.name));
    disp(T)


    %% --- Plot 1D-corrected residuals ---
    figure('Color','w','Position',[100 100 1700 600]);
    b1 = bar(x, y_corr_1D_total, 1); hold on;
    for k = 1:num_datasets, b1(k).FaceColor = colors{k}; end
    set(gca,'XTick',x,'XTickLabel',stn_id,'FontSize',14);
    xtickangle(45);
    xlabel('Station ID'); ylabel('Residual (mm/yr)');
    title('Residuals (1D GIA-corrected total signal): GNSS − Model');
    legend(cellfun(@(d)d.name,datasets,'UniformOutput',false),'Location','northwest', 'FontSize', 14,'Interpreter','none');
    grid on; box on;
    ylim([-7 13]);
    if save_fig
        saveas(gcf, fullfile(fpath_results_figures, sprintf('Residuals_GIA1D_total_corrected_%s_%s.png', rheology_choice, label)));
    end

    %% --- Plot 3D-corrected residuals ---
    figure('Color','w','Position',[100 100 1700 600]);
    b2 = bar(x, y_corr_3D_total, 1); hold on;
    for k = 1:num_datasets, b2(k).FaceColor = colors{k}; end
    set(gca,'XTick',x,'XTickLabel',stn_id,'FontSize',14);
    xtickangle(45);
    xlabel('Station ID'); ylabel('Residual (mm/yr)');
    title('Residuals (3D GIA-corrected total signal): GNSS − Model');
    legend(cellfun(@(d)d.name,datasets,'UniformOutput',false),'Location','northwest', 'FontSize', 14,'Interpreter','none');
    ylim([-7 13]);
    grid on; box on;
    if save_fig
        saveas(gcf, fullfile(fpath_results_figures, sprintf('Residuals_GIA3D_total_corrected_%s_%s.png', rheology_choice, label)));
    end

    %% --- Plot 3D-corrected residuals for Deglacial and Peripheral Glaciers components only---
    figure('Color','w','Position',[100 100 1700 600]);
    b2 = bar(x, y_corr_3D_DG_PG, 1); hold on;
    for k = 1:num_datasets, b2(k).FaceColor = colors{k}; end
    set(gca,'XTick',x,'XTickLabel',stn_id,'FontSize',14);
    xtickangle(45);
    xlabel('Station ID'); ylabel('Residual (mm/yr)');
    title('Residuals (3D GIA-corrected for DG and PG signals): GNSS − Model');
    legend(cellfun(@(d)d.name,datasets,'UniformOutput',false),'Location','northwest', 'FontSize', 14,'Interpreter','none');
    ylim([-7 13]);
    grid on; box on;
    if save_fig
        saveas(gcf, fullfile(fpath_results_figures, sprintf('Residuals_GIA3D_DG_PG_corrected_%s_%s.png', rheology_choice, label)));
    end

    %% --- Plot residuals attributed to LIA---
    figure('Color','w','Position',[100 100 1700 600]);
    b2 = bar(x, gia_val_LIA, 1); hold on;
    set(gca,'XTick',x,'XTickLabel',stn_id,'FontSize',14);
    xtickangle(45);
    xlabel('Station ID'); ylabel('VLM(mm/yr)');
    title('GIA signals that can be attributed to LIA');
    legend('Total 3D GIA minus DG and PG components', 'Location','northwest', 'FontSize', 14,'Interpreter','none');
    grid on; box on;
    ylim([-7 13]);
    if save_fig
        saveas(gcf, fullfile(fpath_results_figures, sprintf('GIA_signals_attributed_to_LIA_%s_%s.png', rheology_choice, label)));
    end

    %% --- Plot residuals corrected for 3D GIA, LIA, and contemporary PG (GrCan). Figure S4b---
    figure('Color','w','Position',[100 100 1700 600]);
    b2 = bar(x, y_corr_total, 1); hold on;
    for k = 1:num_datasets, b2(k).FaceColor = colors{k}; end
    set(gca,'XTick',x,'XTickLabel',stn_id,'FontSize',14);
    xtickangle(45);
    xlabel('Station ID'); ylabel('Residual (mm/yr)');
    title('Residuals (corrected for 3D GIA, LIA, PG and contemporary PG in Canada)');
    legend(cellfun(@(d)d.name,datasets,'UniformOutput',false),'Location','northwest', 'FontSize', 14,'Interpreter','none');
    ylim([-7 13]);
    grid on; box on;
    if save_fig
        saveas(gcf, fullfile(fpath_results_figures, sprintf('Residuals_corrected_for_3D_GIA_LIA_PG_and_contemporary_PG_Can_%s_%s.png', rheology_choice, label)));
    end
end

if any(steps==11)
    % Figure 1
    % Generate a two panel figure of spatial map of mean of total masked ice thickness change (panel 1)
    % and 1-sigma on the ISSM Greenland mesh (panel 2)
    % Note to myself: Covered period varies depending on the availability of ice altimetry product. But for the elastic calculation, it makes sense to show the common period of the GNSS duration (which also varies by stations)

    % spatial map of mean of total masked thickness change on the ISSM mesh
    load_md = true;
    label = 'with_mask';
    plot_stn_id = false;
    plot_total_thickness_change_for_each_ice_data = false;


    if load_md
        mds = cell(numel(datasets), 1);
        for k = 1:numel(datasets)
            ds = datasets{k};
            fname = fullfile(fpath_mesh_model_regional_refined, sprintf('md%d_global_%s.mat', k, label));
            disp(sprintf('loading md data for %s', ds.name));
            mds{k} = loadmodel(fname);
        end
    end

    dice_total = cell(numel(datasets), 1);
    mean_dice_rate = cell(numel(datasets), 1);
    h = cell(numel(datasets), 1);

    % First, choose a common time window for comparison
    common_year_start = 2003;
    common_year_end =2020;
    common_years_total = common_year_end - common_year_start;
    fprintf('Common period selected: %d-%d', common_year_start, common_year_end);


    for k = 1:numel(datasets)
        ds = datasets{k};
        md = mds{k};
        fprintf('\n=== Calculating mean and 1-sigma of total masked ice thickness change in (%s) ===\n', ds.name);
        idx_start = find(md.masstransport.spcthickness(end,:) == common_year_start);
        idx_end = find(md.masstransport.spcthickness(end,:) == common_year_end);
        fprintf('Selecting ice thickness during the common time between year %d and %d',md.masstransport.spcthickness(end,idx_start),md.masstransport.spcthickness(end,idx_end));
        h{k} = md.masstransport.spcthickness(1:end-1,idx_start:idx_end);
        dice = diff(h{k}, 1, 2); % take differences along the time dimension
        dice_total{k} = sum(dice, 2);
        mean_dice_rate{k} = dice_total{k} / common_years_total; % annual trend for each product
    end

    % get mean and sigma across the whole data products
    mean_dice_rate_mat = cat(2, mean_dice_rate{:});  % concatenate along 2nd dimension

    % Compute mean and standard deviation across products (dimension 2)
    mean_dice_rate_all = mean(mean_dice_rate_mat, 2, 'omitnan');
    std_dice_rate_all  = std(mean_dice_rate_mat, 0, 2, 'omitnan');

    % Statistics across the whole ice sheet
    fprintf('Mean of mean_dice_rate_all: %.3f m/yr\n', mean(mean_dice_rate_all, 'omitnan'))
    fprintf('Std of mean_dice_rate_all: %.3f m/yr\n', std(mean_dice_rate_all, 'omitnan'))
    fprintf('Min of mean_dice_rate_all: %.3f m/yr\n', min(mean_dice_rate_all));
    fprintf('Max of mean_dice_rate_all: %.3f m/yr\n', max(mean_dice_rate_all));
    fprintf('Range of mean_dice_rate_all: %.3f m/yr\n', range(mean_dice_rate_all));
    fprintf('Median std_dice_rate_all: %.3f m/yr\n', median(std_dice_rate_all, 'omitnan'))
    fprintf('Min of std_dice_rate_all: %.3f m/yr\n', min(std_dice_rate_all));
    fprintf('Max of std_dice_rate_all: %.3f m/yr\n', max(std_dice_rate_all));
    fprintf('Range of std_dice_rate_all: %.3f m/yr\n', range(range(std_dice_rate_all, 'omitnan')))

    % Panel 1
    disp('Creating Figure 1a')
    % Project to polar stereographic north
    [xg_psn, yg_psn] = ll2psn(lat_gnss, lon_gnss);
    % Extract greenland mesh only
    [meshG, indexG] = extractGreenlandMeshfromGlobal(md);
    data = mean_dice_rate_all(indexG);
    figure('Color','w')
    trisurf(meshG.elements, ...
            meshG.x, meshG.y, ...
            zeros(meshG.numberofvertices,1), ...
            data, ...
            'EdgeColor','none', ...
            'FaceColor','interp');
    view(2); axis equal tight;
    hold on
    caxis([-1 1]);
    greenland('Color',[0.6 0.6 0.6],'LineWidth',0.5);
    colormap(flip(custom_cmap))
    grid off
    axis off
    plot(xg_psn, yg_psn, 'ko', 'MarkerFaceColor','w', 'MarkerSize',6)
    title({'Inter-product mean annual ice thickness change', 'between 2003-2020 (m/yr)'}, 'FontSize',16)
    colorbar()
    %graticulepsn(meshG.lat,meshG.long);
    if save_fig
        exportgraphics(gcf, fullfile(fpath_results_figures,'Figure1a.png'), 'Resolution',300);
    end

    % Panel 2
    disp('Creating Figure 1b')
    data = std_dice_rate_all(indexG);
    figure('Color','w')
    trisurf(meshG.elements, ...
            meshG.x, meshG.y, ...
            zeros(meshG.numberofvertices,1), ...
            data, ...
            'EdgeColor','none', ...
            'FaceColor','interp');
    view(2); axis equal tight;
    hold on
    caxis([0 1]);
    greenland('Color',[0.6 0.6 0.6],'LineWidth',0.5);
    colormap(flip(colormap('hot')))
    grid off
    axis off
    plot(xg_psn, yg_psn, 'ko', 'MarkerFaceColor','w', 'MarkerSize',6)
    title('Inter-product standard deviation (m/yr)', 'FontSize',16)
    colorbar()
    if save_fig
        exportgraphics(gcf, fullfile(fpath_results_figures,'Figure1b.png'), 'Resolution',300);
    end

    % ===== Below plotting is using plotmodel and on a 3D global surface ========
    % panel 1
    %plotmodel(md, 'data', mean_dice_rate_all, 'figure', 10, 'caxis', [-1 1], 'colormap', flip(custom_cmap))
    %hold on
    %title('Inter-product mean thickness change (m/yr)')
    %set(gca,'clipping','off')

    % panel 2
    %xg = r_earth * cosd(lat_gnss) .* cosd(lon_gnss);
    %yg = r_earth * cosd(lat_gnss) .* sind(lon_gnss);
    %zg = r_earth * sind(lat_gnss);

    %plotmodel(md, 'data', std_dice_rate_all, 'figure', 11, 'caxis', [0 0.5], 'colormap',flip(colormap('hot')),'title','1-sigma across datasets')
    %hold on
    %plot3(xg, yg, zg, 'b.', 'MarkerSize', 14);
    %close(1);
    %for i = 1:length(stn_id)
    %    text(xg(i), yg(i), zg(i), stn_id{i}, ...
    %        'FontSize',8, 'Color','b', ...
    %        'VerticalAlignment','bottom');
    %end
    %axis equal off
    %light; lighting gouraud; material dull % add soft 3-D lighting
    %view(53,73)
    %title('Inter-product standard deviation (m/yr)')
    %set(gca,'clipping','off')
    %set(findobj(gca,'Type','text'), 'RotationMode','auto');
    % ===========================================================================

    if plot_total_thickness_change_for_each_ice_data
        % Calculate and plot the maximum and minimum total and mean ice thickness changes for each altimetry data
        % Plot the base data map
        max_total_val = zeros(1,length(datasets));
        idx_total_max = zeros(1,length(datasets));
        min_total_val = zeros(1,length(datasets));
        idx_total_min = zeros(1,length(datasets));
        max_mean_val = zeros(1,length(datasets));
        idx_mean_max = zeros(1,length(datasets));
        min_mean_val = zeros(1,length(datasets));
        idx_mean_min = zeros(1,length(datasets));

        for k = 1:length(datasets)
            ds = datasets{k};
            md = mds{k};

            vals_total = dice_total{k}(indexG);
            vals_mean = mean_dice_rate{k}(indexG);

            % Find the index and value of the maximum
            [max_total_val(k), idx_total_max(k)] = max(vals_total);
            [min_total_val(k), idx_total_min(k)] = min(vals_total);
            [max_mean_val(k), idx_mean_max(k)] = max(vals_mean);
            [min_mean_val(k), idx_mean_min(k)] = min(vals_mean);

            % Extract coordinates of the max and min vertices
            x_max = meshG.x(idx_total_max(k));
            y_max = meshG.y(idx_total_max(k));
            x_min = meshG.x(idx_total_min(k));
            y_min = meshG.y(idx_total_min(k));

            figure('Color','w')
            trisurf(meshG.elements, ...
                    meshG.x, meshG.y, ...
                    zeros(meshG.numberofvertices,1), ...
                    vals_total, ...
                    'EdgeColor','none', ...
                    'FaceColor','interp');
            view(2); axis equal tight;
            hold on
            caxis([-100 100]);
            greenland('Color',[0.6 0.6 0.6],'LineWidth',0.5);
            colormap(flip(custom_cmap))
            grid off
            axis off
            plot(xg_psn, yg_psn, 'ko', 'MarkerFaceColor','w', 'MarkerSize',6)
            plot(x_max, y_max, 'k+', 'MarkerSize', 6, 'LineWidth', 2)
            plot(x_min, y_min, 'k*', 'MarkerSize', 6, 'LineWidth', 2)
            title(sprintf(['Total thickness change between %d and %d \nfor %s (m)'], common_year_start, common_year_end, ds.name))
            colorbar()

            % Optionally add a label
            %text(x_max, y_max, sprintf('Max total: %.2f m \nMax mean: %.2f m/yr', max_total_val(k), max_mean_val(k)), ...
            %     'Color', 'k', 'FontSize', 7, 'FontWeight', 'bold', ...
            %     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
            % Optionally add a label
            %text(x_min, y_min, sprintf('Min total: %.2f m \nMin mean: %.2f m/yr, ', min_total_val(k), min_mean_va(k))), ...
            %     'Color', 'k', 'FontSize', 7, 'FontWeight', 'bold', ...
            %     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
            sprintf('For dataset: %s', ds.name)
            sprintf('Max total change: %.2f m \nMax mean change: %.2f m/yr at index %d' , max_total_val(k), max_mean_val(k), idx_mean_max(k))
            sprintf('Min total change: %.2f m \nMin mean change: %.2f m/yr at index %d' , min_total_val(k), min_mean_val(k), idx_mean_min(k))

            if plot_stn_id
                for i = 1:length(stn_id)
                    text(xg_psn(i), yg_psn(i), stn_id{i}, ...
                        'FontSize',10, 'FontWeight','bold', ...
                        'Color','k', ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','bottom');
                end
            end
        end
    end
end

if any(steps==12)

    % Generates panels for the corresponding supplementary figure (Figure S3)
    % Requires running Step(fig1) in advance.

    save_fig = false;
    plot_stn_id = false;
    diff_product = zeros(meshG.numberofvertices, numel(datasets));
    % Calculate and plot differences bewteen mean and each dataset in total thickness change rate(Supp. Fig. 1)
    for i = 1:length(datasets)
        [maxval, idxmax] = max(mean_dice_rate{i}(indexG) - mean_dice_rate_all(indexG));
        [minval, idxmin] = min(mean_dice_rate{i}(indexG) - mean_dice_rate_all(indexG));

        % Extract coordinates of the max vertex (from your earlier result)
        x_max = meshG.x(idxmax);
        y_max = meshG.y(idxmax);
        x_min = meshG.x(idxmin);
        y_min = meshG.y(idxmin);

        data = mean_dice_rate{i}(indexG) - mean_dice_rate_all(indexG);
        diff_product(:,i) = data;

        figure('Color','w')
        trisurf(meshG.elements, ...
                meshG.x, meshG.y, ...
                zeros(meshG.numberofvertices,1), ...
                data, ...
                'EdgeColor','none', ...
                'FaceColor','interp');
        view(2); axis equal tight;
        hold on
        caxis([-1 1]);
        greenland('Color',[0.6 0.6 0.6],'LineWidth',0.5);
        colormap(flip(custom_cmap))
        grid off
        axis off
        plot(xg_psn, yg_psn, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize',6)
        plot(x_max, y_max, 'k+', 'MarkerSize', 6, 'LineWidth', 2)
        plot(x_min, y_min, 'k*', 'MarkerSize', 6, 'LineWidth', 2)
        title(sprintf(['Deviation from inter-product mean thickness change\nfor %s (m/yr)'], datasets{i}.name), 'FontSize',16,'Interpreter','none')
        colorbar()

        if plot_stn_id
            for i = 1:length(stn_id)
                text(xg_psn(i), yg_psn(i), stn_id{i}, ...
                    'FontSize',10, 'FontWeight','bold', ...
                    'Color','k', ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom');
            end
        end

        % --- Add text showing min and max values in panel ---
        xlims = xlim;
        ylims = ylim;
        text(xlims(1) + 0.02*range(xlims), ylims(2) - 0.05*range(ylims), ...
            sprintf('+ max = %.2f m/yr', maxval), 'Color', 'k', 'FontSize', 15, 'FontWeight', 'bold')

        text(xlims(1) + 0.02*range(xlims), ylims(2) - 0.10*range(ylims), ...
            sprintf('* min = %.2f m/yr', minval), 'Color', 'k', 'FontSize', 15, 'FontWeight', 'bold')

        if save_fig
            saveas(gcf, fullfile(fpath_results_figures, sprintf('Supp_Fig3_panel%d.png', i)));
        end
    end

    % Debug: check if diff_product sums to zero
    diff_sum = sum(diff_product,2); %sum across the product
    sprintf('min_diff_sum: %0.3e and max_diff_sum %0.3e',min(diff_sum), max(diff_sum));
end

if any(steps==13)
    % Generate a figure equivalent to Fig 1 but for VLM
    % Note: This step assumes all the required variables to calculate the mean and std fields have been loaded.
    %       While Fig 1 covers time period (2003-2020) where all altimetry products overlap,
    %       VLM rates shown in this figure are based on different time period for each GNSS station

    save_fig = false;
    plot_stn_id = false;

    % First, do some preprocessing for the modeled elastic VLM rates for the measureItsLive data such that
    % the elastic GrPG signals from Berg et al (2024) is added. This is because measureItsLive data already contains
    % elastic PG signals. Later, we aggregate the modeled elastic VLM with the elastic GrPG signals using Berg et al. (2024),
    % but this leads to double counting the counting elastic GrPG signals in the VLM, which we want to avoid.

    rates_matrix_corr = rates_matrix;
    for i = 1:num_stations
        name_gnss = string(stn_id{i});
        match_idx = find(elas_data.station == name_gnss, 1);
        if ~isempty(match_idx)
            rates_matrix_corr(i,1) = rates_matrix(i,1) - elas_data.Uelastic_GrPG(match_idx);
            rates_matrix_corr(i,2) = rates_matrix(i,2) - elas_data.Uelastic_GrPG(match_idx);
        else
            warning('⚠️ No GIA match found for station %s. Setting corrected residual value to be NaN', name_gnss);
        end
    end

    % For raw residuals between GNSS and model elastic VLM
    mean_vlm_rates = mean(rates_matrix, 2, 'omitnan');       % → T×1 mean across datasets
    std_vlm_rates  = std(rates_matrix,0,2);      % → T×1 sigma across datasets
    vlm_rates_q25 = prctile(rates_matrix, 25, 2);
    vlm_rates_q75 = prctile(rates_matrix, 75, 2);
    vlm_iqr = vlm_rates_q75 - vlm_rates_q25;

    % For redisuals where modeled VLM is corrected for GIA and PGs signals
    mean_vlm_rates_corr = mean(rates_matrix_corr, 2, 'omitnan');       % → T×1 mean across datasets where JPL data are corrected for GrPG
    std_vlm_rates_corr  = std(rates_matrix_corr,0,2);      % → T×1 sigma across datasets
    vlm_rates_q25_corr = prctile(rates_matrix_corr, 25, 2);
    vlm_rates_q75_corr = prctile(rates_matrix_corr, 75, 2);
    vlm_iqr_corr = vlm_rates_q75_corr - vlm_rates_q25_corr;

    T = table(stn_id(:), mean_vlm_rates, std_vlm_rates, ...
    'VariableNames', {'Station','MeanVLM','StdVLM'});

    % panel 1
    disp('Creating Figure 2a')
    figure('Color','w');
    fig = gcf;
    fig.Position(4) = fig.Position(4) + 200;   % +200 pixels vertically

    % --- Base velocity axes ---
    ax1 = axes();
    itslive_imagesc(5);
    hold(ax1,'on');
    greenland('Color',[0.6 0.6 0.6],'LineWidth',0.5);

    ax1.ColorScale = 'log';
    clim(ax1, [1 1e4]);
    view(ax1, 2); axis(ax1,'equal','tight');
    grid(ax1,'off');
    ax1.Visible = 'off';
    colormap(ax1, flip(colormap('hot')));

    % GNSS overlay
    ax2 = axes('Position', ax1.Position, 'Color','none');
    hold(ax2,'on');
    scatter_handle = scatter(ax2, xg_psn, yg_psn, 100, mean_vlm_rates, ...
                         'filled','MarkerEdgeColor','k','LineWidth',0.5);

    colormap(ax2, flip(pink));
    caxis(ax2, [0 15]);

    % Match geometry
    ax2.XLim = ax1.XLim;
    ax2.YLim = ax1.YLim;
    ax2.DataAspectRatio = ax1.DataAspectRatio;
    ax2.PlotBoxAspectRatio = ax1.PlotBoxAspectRatio;
    ax2.XDir = ax1.XDir;
    ax2.YDir = ax1.YDir;
    axis(ax2,'equal','tight');
    ax2.Visible = 'off';

    linkaxes([ax1 ax2],'xy');

    % **Bring GNSS to front**
    uistack(scatter_handle,'top');

    %% --- lock axis positions ---
    ax1.PositionConstraint = 'innerposition';
    ax2.PositionConstraint = 'innerposition';

    %% --- Save original axes position BEFORE adding colorbars
    origPos = ax1.Position;

    %% Add both colorbars
    cb1 = colorbar(ax1,'Location','eastoutside');
    cb1.Label.String = 'Velocity (m/yr, log scale)';

    cb2 = colorbar(ax2,'Location','southoutside');
    cb2.Label.String = 'Mean VLM (mm/yr)';

    %% --- Restore axes positions (critical!)
    ax1.Position = origPos;
    ax2.Position = origPos;

    %% --- Now place colorbars manually
    cb1.Position = [0.8 0.20 0.03 0.60];
    cb2.Position = [0.25 0.077 0.50 0.025];

    % Super-group title
    sgtitle('Mean ice velocity in Greenland and mean VLM at GNSS sites', 'FontSize',18)

    if plot_stn_id
        text(xg_psn, yg_psn + 10e3, stn_id, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',9, 'Color','k');
    end

    if save_fig
        saveas(gcf, fullfile(fpath_results_figures, 'Figure2a.png'));
    end

    % panel 2
    disp('Creating Figure 2b')
    figure('Color','w','Position',[100 100 1500 500]);
    x = 1:num_stations;

    % Add GNSS data points
    % Plot GNSS rates as reference dots
    errorbar(x, gnss_rates_berg, gnss_rates_err_berg, ...
         'o', 'Color','b','MarkerFaceColor','none','MarkerSize',5.5, 'LineWidth',1.2);
    %plot(1:num_stations, gnss_rates, 'ko', 'MarkerSize', 5.5, 'MarkerFaceColor', 'black', 'DisplayName', 'GNSS Rate');
    hold on

    % --- boxchart ---
    Xpos = repmat((1:length(stn_id))', length(datasets), 1);
    Yval = rates_matrix(:);

    % Box-and-whisker style
    bc = boxchart(Xpos, Yval, 'BoxFaceColor', 'k', 'BoxEdgeColor', 'k');
    bc.BoxWidth = 0.6;
    bc.MarkerStyle = 'none';

    % Plot mean VLM and 1-sigma std dev.
    plot(x, mean_vlm_rates, ...
         'o', 'Color','k','MarkerFaceColor','none','MarkerSize',6, 'LineWidth',1.2);
    %errorbar(x, mean_vlm_rates, std_vlm_rates, ...
    %     'o', 'Color','k','MarkerFaceColor','none','MarkerSize',6, 'LineWidth',1.2);
    grid on;
    set(gca,'GridAlpha',0.25);   % lighter grid
    set(gca,'XTick',x,'XTickLabel',stn_id,'FontSize',14);xtickangle(45);
    xlim([0.5 num_stations+0.5]);
    %ylim([0 max(mean_vlm_rates + std_vlm_rates)*1.10]);

    set(gca,'XTick',x,'XTickLabel',stn_id,'FontSize',12); xtickangle(45);
    xlabel('Station ID')
    ylabel('Vertical land motion rate (mm/yr)');
    title('GNSS-observed and modeled elastic VLM rates in at Greenland GNSS stations','FontSize',16)
    legend({'Observed GNSS rate (Berg et al., 2024)','Modeled elastic VLM distribution (median, interquartile range, full range)','Model mean'}, 'Location','northwest');

    if save_fig
        exportgraphics(gcf, fullfile(fpath_results_figures,'Figure2b.png'), 'Resolution',300);
    end

    % Plots for GNSS-Model VLM residuals
    disp('Creating Figure 3a')
    residual = NaN(num_stations, 1);
    residual_err = NaN(num_stations, 1);
    residual = gnss_rates_berg(:) - mean_vlm_rates(:);
    vlm_err = 0.5 * vlm_iqr; % set modeled elastic vlm uncertainty is half of the interquartile range
    residual_err = sqrt(gnss_rates_err_berg(:).^2 + vlm_err(:).^2);
    %residual_err = sqrt(gnss_rates_err_berg.^2 + std_vlm_rates.^2); % this is appropriate when model uncertainty is defined by std dev.

    figure('Color','w','Position',[100 100 1500 400]);
    x = 1:num_stations;
        errorbar(x, residual, residual_err, ...
         'o', 'Color','k','MarkerFaceColor','none','MarkerSize',6, 'LineWidth',1.2);
    hold on
    grid on;
    set(gca,'GridAlpha',0.25);   % lighter grid
    set(gca,'XTick',x,'XTickLabel',stn_id,'FontSize',14);xtickangle(45);
    ylabel('Mean residual and uncertainty (mm/yr)');
    xlabel('Station ID')
    title('Misfits between GNSS-observed and modeled elastic VLM at Greenland GNSS Sites')
    xlim([0.5 53.5]);
    ylim([-6 14])
    yline(0,'k-','LineWidth',1);

    if save_fig
        exportgraphics(gcf, fullfile(fpath_results_figures,'Figure3a.png'), 'Resolution',300);
    end

    % Check
    fprintf('Maximum VLM spread is at station %s with range %.3f mm/yr\n', stn_id{find(residual_err==(max(residual_err)))}, 2*max(residual_err));
    fprintf('Minimum VLM spread is at station %s with range %.3f mm/yr\n', stn_id{find(residual_err==(min(residual_err)))}, 2*min(residual_err));

    % --- Create table of residuals ---
    % Build clean table
    T_residuals = table( ...
        string(stn_id(:)), ...
        round(residual(:), 1), ...
        round(residual_err(:), 1), ...
        'VariableNames', {'Station','Residual (mm/yr)','Sigma(mm/yr)'} ...
    );

    % Display the clean table
    disp(T_residuals)

    % Optional: save
    writetable(T_residuals, fullfile(fpath_results_general,'VLM_residuals_table_raw.csv'));

    disp('Creating Figure 3b') % plot isolating residual only due to GIA, LIA and 20th-Century PGs.
    residual_corrected = NaN(num_stations, 1);
    residual_corrected_err = NaN(num_stations, 1);
    for i = 1:num_stations
        name_gnss = string(stn_id{i});
        match_idx_elas = find(elas_data.station == name_gnss, 1);
        if ~isempty(match_idx_elas) % if the station is in the elas_data table
            residual_corrected(i) = gnss_rates_berg(i) - (mean_vlm_rates_corr(i) + elas_data.Uelastic_CanPG(match_idx_elas) + elas_data.Uelastic_GrPG(match_idx_elas));
            residual_corrected_err(i) = sqrt(gnss_rates_err_berg(i).^2 + (0.5 * vlm_iqr_corr(i)).^2 + elas_data.Uelastic_CanPG_sigma(match_idx_elas).^2+ elas_data.Uelastic_GrPG_sigma(match_idx_elas).^2);
            %residual_corrected_err(i) = sqrt(gnss_rates_err_berg(i).^2 + std_vlm_rates_corr(i).^2 + elas_data.Uelastic_CanPG_sigma(match_idx_elas).^2+ elas_data.Uelastic_GrPG_sigma(match_idx_elas).^2);
        else
            warning('⚠️ No elastic VLM match found for station %s. Setting corrected residual value to be NaN', name_gnss);
        end
    end

    figure('Color','w','Position',[100 100 1500 400]);
    x = 1:num_stations;
        errorbar(x, residual_corrected, residual_corrected_err, ...
         'o', 'Color','k','MarkerFaceColor','none','MarkerSize',6, 'LineWidth',1.2);
    hold on
    grid on;
    set(gca,'GridAlpha',0.25);   % lighter grid
    set(gca,'XTick',x,'XTickLabel',stn_id,'FontSize',14);xtickangle(45);
    ylabel('Mean residual uncertainty (mm/yr)');
    xlabel('Station ID')
    title('Misfits between GNSS-observed and modeled VLM corrected for contemporary PGs signals at Greenland GNSS Sites')
    xlim([0.5 53.5]);
    ylim([-7 14])
    yline(0,'k-','LineWidth',1);

    if save_fig
        exportgraphics(gcf, fullfile(fpath_results_figures,'Figure3b.png'), 'Resolution',300);
    end

    % Build clean table
    T_residuals = table( ...
        string(stn_id(:)), ...
        round(residual_corrected(:), 1), ...
        round(residual_corrected_err(:), 1), ...
        'VariableNames', {'Station','Residual (mm/yr)','Sigma(mm/yr)'} ...
    );

    % Optional: save
    writetable(T_residuals, fullfile(fpath_results_general,'VLM_residuals_table_corrected_for_GrCanPGs.csv'));

    % Check
    fprintf('Maximum VLM spread is at station %s with range %.2f mm/yr\n', stn_id{find(residual_corrected_err==(max(residual_corrected_err)))}, 2*max(residual_corrected_err));
    fprintf('Minimum VLM spread is at station %s with range %.2f mm/yr\n', stn_id{find(residual_corrected_err==(min(residual_corrected_err)))}, 2*min(residual_corrected_err));

    disp('Creating Figure 3c')
    residual_corrected = NaN(num_stations, 1);
    residual_corrected_err = NaN(num_stations, 1);
    for i = 1:num_stations
        name_gnss = string(stn_id{i});
        match_idx_gia = find(gia_data.name == name_gnss, 1);
        match_idx_elas = find(elas_data.station == name_gnss, 1);
        if ~isempty(match_idx_gia) && ~isempty(match_idx_elas)
            residual_corrected(i) = gnss_rates_berg(i) - (mean_vlm_rates_corr(i) + gia3D_total_mean(match_idx_gia) + elas_data.Uelastic_CanPG(match_idx_elas) + elas_data.Uelastic_GrPG(match_idx_elas));
            residual_corrected_err(i) = sqrt(gnss_rates_err_berg(i).^2 + (0.5 * vlm_iqr_corr(i)).^2 + gia3D_total_mean_sigma(match_idx_gia).^2 + elas_data.Uelastic_CanPG_sigma(match_idx_elas).^2+ elas_data.Uelastic_GrPG_sigma(match_idx_elas).^2);
            %residual_corrected_err(i) = sqrt(gnss_rates_err_berg(i).^2 + std_vlm_rates_corr(i).^2 + gia3D_total_mean_sigma(match_idx_gia).^2 + elas_data.Uelastic_CanPG_sigma(match_idx_elas).^2+ elas_data.Uelastic_GrPG_sigma(match_idx_elas).^2);
        else
            warning('⚠️ No GIA match found for station %s. Setting corrected residual value to be NaN', name_gnss);
        end
    end

    figure('Color','w','Position',[100 100 1500 400]);
    x = 1:num_stations;
        errorbar(x, residual_corrected, residual_corrected_err, ...
         'o', 'Color','k','MarkerFaceColor','none','MarkerSize',6, 'LineWidth',1.2);
    hold on
    grid on;
    set(gca,'GridAlpha',0.25);   % lighter grid
    set(gca,'XTick',x,'XTickLabel',stn_id,'FontSize',14);xtickangle(45);
    ylabel('Mean residual uncertainty (mm/yr)');
    xlabel('Station ID')
    title('Misfits between GNSS-observed and modeled total VLM at Greenland GNSS Sites')
    xlim([0.5 53.5]);
    ylim([-6 14])
    yline(0,'k-','LineWidth',1);

    if save_fig
        exportgraphics(gcf, fullfile(fpath_results_figures,'Figure3c.png'), 'Resolution',300);
    end
    % Check
    fprintf('Maximum VLM spread is at station %s with range %.2f mm/yr\n', stn_id{find(residual_corrected_err==(max(residual_corrected_err)))}, 2*max(residual_corrected_err));
    fprintf('Minimum VLM spread is at station %s with range %.2f mm/yr\n', stn_id{find(residual_corrected_err==(min(residual_corrected_err)))}, 2*min(residual_corrected_err));


    % Build clean table
    T_residuals = table( ...
        string(stn_id(:)), ...
        round(residual_corrected(:), 1), ...
        round(residual_corrected_err(:), 1), ...
        'VariableNames', {'Station','Residual (mm/yr)','Sigma(mm/yr)'} ...
    );

    % Optional: save
    writetable(T_residuals, fullfile(fpath_results_general,'VLM_residuals_table_corrected_for_LIA_GIA_AllPGs.csv'));

end

if any(steps==14)
    save_fig = false;
    % Create Figure 4 showing the envolup of residual across Greenland
    % Requires having Steps 10, 11, 12 already run. (need the variable 'y' for raw residual, 'mean_dice_rate' and 'meshG' info)
    for i = 1:length(stn_id)
        [max_residual(i), idx_max(i)] = max(y(i,:), [], 2);
        [min_residual(i), idx_min(i)] = min(y(i,:), [], 2);
        envelope_residual(i) = max_residual(i) - min_residual(i);
    end

    disp('Plotting residual envelope at Greenland GNSS sites') % Not used in the manuscript
    figure('Color','w');
    fig = gcf;
    fig.Position(4) = fig.Position(4) + 200;   % +200 pixels vertically
    itslive_imagesc(5)
    hold on
    greenland('Color',[0.6 0.6 0.6],'LineWidth',0.5);
    set(gca,'colorscale','log')
    clim([1 10e3])
    view(2); axis equal tight;
    greenland('Color',[0.6 0.6 0.6],'LineWidth',0.5);
    colormap(flip(custom_cmap))
    grid off
    ax1 = gca;
    ax1.Visible = 'off';

    % GNSS overlay
    ax2 = axes('Position', ax1.Position, 'Color','none');
    hold(ax2,'on');
    scatter_handle = scatter(ax2, xg_psn, yg_psn, 100, envelope_residual, ...
                         'filled','MarkerEdgeColor','k','LineWidth',0.5);

    colormap(ax2, flip(pink));
    caxis(ax2, [0 8]);

    % Match geometry
    ax2.XLim = ax1.XLim;
    ax2.YLim = ax1.YLim;
    ax2.DataAspectRatio = ax1.DataAspectRatio;
    ax2.PlotBoxAspectRatio = ax1.PlotBoxAspectRatio;
    ax2.XDir = ax1.XDir;
    ax2.YDir = ax1.YDir;
    axis(ax2,'equal','tight');
    ax2.Visible = 'off';

    linkaxes([ax1 ax2],'xy');

    % **Bring GNSS to front**
    uistack(scatter_handle,'top');

    %% --- lock axis positions ---
    ax1.PositionConstraint = 'innerposition';
    ax2.PositionConstraint = 'innerposition';

    %% --- Save original axes position BEFORE adding colorbars
    origPos = ax1.Position;

    %% Add both colorbars
    cb1 = colorbar(ax1,'Location','eastoutside');
    cb1.Label.String = 'Velocity (m/yr, log scale)';

    cb2 = colorbar(ax2,'Location','southoutside');
    cb2.Label.String = 'Residual envelope (mm/yr)';

    %% --- Restore axes positions (critical!)
    ax1.Position = origPos;
    ax2.Position = origPos;

    %% --- Now place colorbars manually
    cb1.Position = [0.8 0.20 0.03 0.60];
    cb2.Position = [0.25 0.077 0.50 0.025];

    % Super-group title
    sgtitle('Residual envelope at GNSS sites', 'FontSize',16)

    if plot_stn_id
        text(xg_psn, yg_psn + 10e3, stn_id, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',9, 'Color','k');
    end

    if save_fig
        exportgraphics(gcf, fullfile(fpath_results_figures,'Residual_envelope_at_Greenland_GNSS_sites.png'), 'Resolution',300);
    end

    disp('Creating Figure 4a')
    figure('Color','w','Position',[100 100 1700 350]);
    b2 = bar(x, envelope_residual, 1); hold on;
    b2.FaceColor = 'blue';
    set(gca,'XTick',x,'XTickLabel',stn_id,'FontSize',16);
    xtickangle(45);
    xlabel('Station ID','FontSize',20); ylabel('Residual envelope (mm/yr)','FontSize',20);
    title('Inter-product residual envelope (max minus min) at Greenland GNSS stations','FontSize',24);
    yline(1,'k--','LineWidth',1);
    grid on; box on;
    if save_fig
        exportgraphics(gcf, fullfile(fpath_results_figures,'Figure4a.png'), 'Resolution',300);
    end

    % Plot mean ice elevation change rate difference between limiting altimetry-FAC models
    disp('Creating Figure 4')
    stn_names = {'HEL2','UTMG','DKSG'};
    pad_km = 100;                 % zoom half-width in km
    pad_m  = pad_km * 1000;
    for i=1:length(stn_id)
        if any(strcmp(stn_id{i}, stn_names))
            % find which ice models are limiting models by looking up residuals
            fprintf('%s: model that gives min residual =%s, max residual=%s\n', ...
                stn_id{i}, dataset_names{idx_min(i)}, dataset_names{idx_max(i)}); % dataset that gives the min and max residuals for the given station
            data = mean_dice_rate{idx_min(i)}(indexG) - mean_dice_rate{idx_max(i)}(indexG);

            figure('Color','w')
            trisurf(meshG.elements, ...
                    meshG.x, meshG.y, ...
                    zeros(meshG.numberofvertices,1), ...
                    data, ...
                    'EdgeColor','none', ...
                    'FaceColor','interp');
            view(2); axis equal tight;
            hold on
            clim([-10 10]);
            greenland('Color',[0.6 0.6 0.6],'LineWidth',0.7);
            colormap(flip(custom_cmap))
            grid off
            axis off
            plot(xg_psn, yg_psn, 'ko', 'MarkerFaceColor','w', 'MarkerSize',8)

            % ---- ZOOM around the selected station ----
            x0 = xg_psn(i);
            y0 = yg_psn(i);
            xlim([x0 - pad_m, x0 + pad_m]);
            ylim([y0 - pad_m, y0 + pad_m]);

            title(sprintf(['Difference in mean ice elevation change (m/yr)\n' ...
                'at %s between %s and %s'], ...
                stn_id{i}, datasets{idx_min(i)}.name, datasets{idx_max(i)}.name), ...
                'FontSize',16)

            colorbar()

            % label ONLY the stations that fall inside the zoom window
            inwin = xg_psn >= x0-pad_m & xg_psn <= x0+pad_m & ...
                    yg_psn >= y0-pad_m & yg_psn <= y0+pad_m;

            for k = find(inwin(:))'
                text(xg_psn(k), yg_psn(k), stn_id{k}, ...
                    'FontSize',13, 'FontWeight','bold', ...
                    'Color','k', 'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom');
            end

            if save_fig
                exportgraphics(gcf, fullfile(fpath_results_figures,sprintf('Figure4_%s_diff_in_mean_ice_elevation_change_2003_2020.png', stn_id{i})), 'Resolution',300);
            end

            % ====== Plot mean ice elevation change in the same spatial domain ======

            data = mean_dice_rate_all(indexG) ;

            figure('Color','w')
            trisurf(meshG.elements, ...
                    meshG.x, meshG.y, ...
                    zeros(meshG.numberofvertices,1), ...
                    data, ...
                    'EdgeColor','none', ...
                    'FaceColor','interp');
            view(2); axis equal tight;
            hold on
            clim([-10 10]);
            greenland('Color',[0.6 0.6 0.6],'LineWidth',0.7);
            colormap(flip(custom_cmap))
            grid off
            axis off
            plot(xg_psn, yg_psn, 'ko', 'MarkerFaceColor','w', 'MarkerSize',8)
            title('Mean ice elevation change (m/yr) across all products')
            % ---- ZOOM around the selected station ----
            x0 = xg_psn(i);
            y0 = yg_psn(i);
            xlim([x0 - pad_m, x0 + pad_m]);
            ylim([y0 - pad_m, y0 + pad_m]);

            colorbar()

            % label ONLY the stations that fall inside the zoom window
            inwin = xg_psn >= x0-pad_m & xg_psn <= x0+pad_m & ...
                    yg_psn >= y0-pad_m & yg_psn <= y0+pad_m;

            for k = find(inwin(:))'
                text(xg_psn(k), yg_psn(k), stn_id{k}, ...
                    'FontSize',13, 'FontWeight','bold', ...
                    'Color','k', 'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom');
            end

            if save_fig
                exportgraphics(gcf, fullfile(fpath_results_figures,sprintf('Figure4_%s_interproduct_mean_ice_elevation_change.png', stn_id{i})), 'Resolution',300);
            end
        end
    end
end
