function [ice_masks, years, lat_mask, lon_mask, x_mask, y_mask] = preprocess_glacier_mask(target_years)
    % preprocess_glacier_mask.m
    % Holly Han (created: July 29th, 2025; Last edited: August 19th, 2025).
    % Processes Chad Green's glacier mask datasets.
    % Data source: https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/MEASURES/NSIDC-0793/1/1972/09/15/NSIDC-0793_19720915-20220215_V01.0.nc
    %
    % Inputs:
    %   - target_years: Target years to extract (optional, default: 1993-2022)
    % Outputs:
    %    - ice_masks: Ice mask data [lat x lon x years]
    %    - years: Year vector
    %    - lat_mask: Latitude coordinates on sphere
    %    - lon_mask: Longitude coordinates on sphere
    %    - x_mask: X coordinates in projected system (m)
    %    - y_mask: Y coordinates in projected system (m)
    
    % Set default parameters
    if nargin < 1
        target_years = 1993:2022; % Default years to extract
    end
    
    % File path to Chad Green's glacier mask dataset
    mask_file = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/glacier_mask/greenland_ice_masks_1972-2022_v1.nc';
    
    % Check if file exists
    if ~exist(mask_file, 'file')
        error('Glacier mask file not found: %s', mask_file);
    end
    
    % Read the dataset
    try
        % Read coordinates
        x_orig = ncread(mask_file, 'x');  % X coordinates in meters (EPSG:3413)
        y_orig = ncread(mask_file, 'y');  % Y coordinates in meters (EPSG:3413)
        time_data = ncread(mask_file, 'time'); % Time variable
        
        % Handle time conversion - time is "days since 1900-1-1 0:0:0"
        % Convert days since 1900-01-01 to years
        ref_date = datetime(1900,1,1);
        dates_all = ref_date + days(time_data);
        years_all = year(dates_all);
        doy_all = day(dates_all,'dayofyear');
        target_doy = 197; % ~July 15
        
        time_indices = [];
        years = [];
        
        for yr = target_years
            idxs = find(years_all == yr);
            if isempty(idxs), continue, end
            [~, mid_idx] = min(abs(doy_all(idxs) - target_doy));
            time_indices(end+1) = idxs(mid_idx);
            years(end+1) = yr;
        end
        
        disp(['Selected ', num2str(length(years)), ' time slices from ', num2str(years(1)), ' to ', num2str(years(end))]);
        disp(['Time indices to read: ', num2str(time_indices)]);
        
        % Check if our indices are valid
        if max(time_indices) > length(time_data)
            error('Time indices exceed file dimensions! Max index: %d, File size: %d', max(time_indices), length(time_data));
        end
        
        % Read ice mask data for selected time slices only
        if length(time_indices) > 10
            disp('Reading ice mask data in chunks...');
            ice_data = [];
            chunk_size = 5;
            
            for i = 1:chunk_size:length(time_indices)
                end_idx = min(i + chunk_size - 1, length(time_indices));
                chunk_indices = time_indices(i:end_idx);
                disp(['Reading chunk ', num2str(i), '-', num2str(end_idx), ' (years ', num2str(years(i)), '-', num2str(years(end_idx)), ')']);
                start_vals = [1, 1, chunk_indices(1)];
                count_vals = [Inf, Inf, length(chunk_indices)];
                chunk_data = ncread(mask_file, 'ice', start_vals, count_vals);
                if isempty(ice_data)
                    ice_data = chunk_data;
                else
                    ice_data = cat(3, ice_data, chunk_data);
                end
            end
        else
            disp('Reading selected ice mask data individually...');
            ice_data = [];
            info = ncinfo(mask_file, 'ice');
            for i = 1:length(time_indices)
                disp(['Reading time slice ', num2str(i), ' (index ', num2str(time_indices(i)), ', year ', num2str(years(i)), ')']);
                start_vals = [1, 1, time_indices(i)];
                count_vals = [Inf, Inf, 1];
                single_slice = ncread(mask_file, 'ice', start_vals, count_vals);
                disp(['Single slice size: ', num2str(size(single_slice))]);
                if isempty(ice_data)
                    ice_data = single_slice;
                else
                    ice_data = cat(3, ice_data, single_slice);
                end
            end
        end
        
        disp(['Loaded mask data from ', num2str(years(1)), ' to ', num2str(years(end))]);
        disp(['ice_data size (raw): ', num2str(size(ice_data))]);
        disp(['Years requested: ', num2str(years)]);
        
    catch ME
        disp('Available variables in the NetCDF file:');
        info = ncinfo(mask_file);
        for i = 1:length(info.Variables)
            disp(['  - ', info.Variables(i).Name]);
        end
        rethrow(ME);
    end
    
    [X_2d, Y_2d] = meshgrid(x_orig, y_orig);
    x_mask = X_2d;
    y_mask = Y_2d;
    
    disp('Converting coordinates from EPSG:3413 to geographic...');
    proj_info = projcrs(3413);
    [lat_mask, lon_mask] = projinv(proj_info, X_2d, Y_2d);
    
    ice_masks = permute(ice_data, [2, 1, 3]);
    disp(['ice_masks size after permute: ', num2str(size(ice_masks))]);
    
    ice_masks = ice_masks > 0;
    ice_masks(isnan(ice_masks)) = false;
    ice_masks = double(ice_masks);

    disp(['Processed mask data for years: ', num2str(years(1)), '-', num2str(years(end))]);
    disp(['Final mask size: ', num2str(size(ice_masks))]);
    
    disp('Glacier mask preprocessing completed!');
end
