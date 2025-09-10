function [ice_masks, years, lat_mask, lon_mask, x_mask, y_mask] = preprocess_glacier_mask(target_years)
    % preprocess_glacier_mask.m
    % Holly Han (created: July 29th, 2025; Last edited: August 19th, 2025).
    % Processes Chad Green's glacier mask datasets.
    % Data source: https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/MEASURES/NSIDC-0793/1/1972/09/15/NSIDC-0793_19720915-20220215_V01.0.nc
    %
    % Inputs:
    %   - target_years: Target years to extract (optional, default: 1995-1996)
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
        days_since_1900 = time_data;
        years_since_1900 = days_since_1900 / 365.25; % Account for leap years
        years_all = floor(1900 + years_since_1900);
        
        time_indices = [];
        years = [];
        
        for year = target_years
            year_idx = find(years_all == year, 1, 'first'); % Get first occurrence of each year
            if ~isempty(year_idx)
                time_indices(end+1) = year_idx;
                years(end+1) = year;
            end
        end
        
        disp(['Selected ', num2str(length(years)), ' time slices from ', num2str(years(1)), ' to ', num2str(years(end))]);
        disp(['Time indices to read: ', num2str(time_indices)]);
        
        % Debug time information
        %disp(['Total time values in file: ', num2str(length(days_since_1900))]);
        %disp(['Year range in file: ', num2str(min(years_all)), ' to ', num2str(max(years_all))]);
        %disp(['Max time index we want to read: ', num2str(max(time_indices))]);
        
        % Check if our indices are valid
        if max(time_indices) > length(days_since_1900)
            error('Time indices exceed file dimensions! Max index: %d, File size: %d', max(time_indices), length(days_since_1900));
        end
        
        % Read ice mask data for selected time slices only
        if length(time_indices) > 10
            % If still too many slices, read in chunks
            disp('Reading ice mask data in chunks...');
            ice_data = [];
            chunk_size = 5; % Process 5 years at a time
            
            for i = 1:chunk_size:length(time_indices)
                end_idx = min(i + chunk_size - 1, length(time_indices));
                chunk_indices = time_indices(i:end_idx);
                
                disp(['Reading chunk ', num2str(i), '-', num2str(end_idx), ' (years ', num2str(years(i)), '-', num2str(years(end_idx)), ')']);
                
                % Read this chunk - ice variable is ice(x, y, time)
                start_vals = [1, 1, chunk_indices(1)];
                count_vals = [Inf, Inf, length(chunk_indices)];
                chunk_data = ncread(mask_file, 'ice', start_vals, count_vals);
                
                if isempty(ice_data)
                    ice_data = chunk_data;
                else
                    ice_data = cat(3, ice_data, chunk_data); % Concatenate along time dimension (3rd)
                end
            end
        else
            % Read selected time slices individually - ice variable is ice(time, y, x)
            disp('Reading selected ice mask data individually...');
            ice_data = [];
            
            % First, let's check the actual dimensions of the ice variable
            info = ncinfo(mask_file, 'ice');
            %disp(['Ice variable dimensions: ', num2str(info.Size)]);
            %disp(['Expected: [time=594, y=23334, x=13333]']);
            %disp(['Coordinate lengths: x=', num2str(length(x_orig)), ', y=', num2str(length(y_orig))]);
            
            for i = 1:length(time_indices)
                disp(['Reading time slice ', num2str(i), ' (index ', num2str(time_indices(i)), ', year ', num2str(years(i)), ')']);
                
                % The ice variable is [x, y, time] = [13333, 23334, 594]
                % So to read one time slice, we need: start=[1, 1, time_index], count=[Inf, Inf, 1]
                start_vals = [1, 1, time_indices(i)];
                count_vals = [Inf, Inf, 1];
                single_slice = ncread(mask_file, 'ice', start_vals, count_vals);
                disp(['Single slice size: ', num2str(size(single_slice))]);
                
                if isempty(ice_data)
                    ice_data = single_slice;
                else
                    ice_data = cat(3, ice_data, single_slice); % Concatenate along time dimension (3rd)
                end
            end
        end
        
        disp(['Loaded mask data from ', num2str(years(1)), ' to ', num2str(years(end))]);
        disp(['ice_data size (raw): ', num2str(size(ice_data))]);
        
        % Check if we have the right coordinate dimensions
        %disp(['Expected dimensions based on coordinates:']);
        %disp(['  x coordinates: ', num2str(length(x_orig))]);
        %disp(['  y coordinates: ', num2str(length(y_orig))]);
        %disp(['  time indices selected: ', num2str(length(time_indices))]);
        %disp(['  Expected ice data: [time, y, x] = [', num2str(length(time_indices)), ', ', num2str(length(y_orig)), ', ', num2str(length(x_orig)), ']']);
        
        % Check what we actually read
        %disp(['Actual ice_data size: ', num2str(size(ice_data))]);
        
        % The issue might be that we have too many time indices
        %disp(['Time indices used: ', num2str(time_indices)]);
        disp(['Years requested: ', num2str(years)]);
        
    catch ME
        % Try alternative variable names or provide more specific error
        disp('Available variables in the NetCDF file:');
        info = ncinfo(mask_file);
        for i = 1:length(info.Variables)
            disp(['  - ', info.Variables(i).Name]);
        end
        rethrow(ME);
    end
    
    % Create 2D coordinate arrays (meshgrid)
    % x_orig has 13333 elements, y_orig has 23334 elements
    % meshgrid creates [length(y_orig), length(x_orig)] = [23334, 13333]
    [X_2d, Y_2d] = meshgrid(x_orig, y_orig);
    
    % But ice_masks is [23334, 594, 3] which doesn't match
    % This suggests we need to check what's happening with the ice data reading
    
    % Store original projected coordinates
    x_mask = X_2d;
    y_mask = Y_2d;
    
    % Convert to geographic coordinates (EPSG:3413 to WGS84)
    disp('Converting coordinates from EPSG:3413 to geographic...');
    proj_info = projcrs(3413); % Polar Stereographic North
    [lat_ellipsoid, lon_ellipsoid] = projinv(proj_info, X_2d, Y_2d);
    
    % Transform to sphere for consistency with other datasets
    r_earth = 6371000; % Earth radius in meters
    [lat_mask, lon_mask, ~] = ellipsoid_to_sphere(lat_ellipsoid, lon_ellipsoid, r_earth, zeros(size(lat_ellipsoid)), false);
    
    % Process ice mask data
    % ice_data is ice(time, y, x), we want [y, x, time] for MATLAB convention
    ice_masks = permute(ice_data, [2, 1, 3]);
    
    % Debug: check dimensions
    disp(['ice_masks size after permute: ', num2str(size(ice_masks))]);
    %disp(['x_orig size: ', num2str(size(x_orig))]);
    %disp(['y_orig size: ', num2str(size(y_orig))]);
    %disp(['X_2d size: ', num2str(size(X_2d))]);
    %disp(['Y_2d size: ', num2str(size(Y_2d))]);
    
    % Data already filtered to target years during reading
    
    % Apply basic processing
    % Convert to logical mask (ice = 1, no ice = 0)
    ice_masks = ice_masks > 0;
    
    % Set no-data values
    ice_masks(isnan(ice_masks)) = false;
    
    % Convert to double
    ice_masks = double(ice_masks);

    disp(['Processed mask data for years: ', num2str(years(1)), '-', num2str(years(end))]);
    disp(['Final mask size: ', num2str(size(ice_masks))]);
    
    disp('Glacier mask preprocessing completed!');
end 