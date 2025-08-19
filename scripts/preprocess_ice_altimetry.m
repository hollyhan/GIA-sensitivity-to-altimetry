function [h_annual, dhdt_annual, dhdt_monthly, years_thickness, lat_sphere, long_sphere, X, Y] = preprocess_ice_altimetry(data_name, plot_altimetry)
    % preprocess_ice_altimetry.m
    % Holly Han (created: July 25th, 2025; Last edited: August 13th, 2025).
    % Preprocesses ice thickness and elevation change data from altimetry.
    %
    % Inputs:
    %   - data_name: Dataset name ('measureItsLive', 'DTU2016', 'DTU2025', 'Buffalo2025-GEMB', 'Buffalo2025-GSFC', 'Buffalo2025-IMAU')
    %   - plot_altimetry: boolean to plot the altimetry data (true or false)
    %
    % Outputs:
    %    - h_annual: ice thickness at annual time interval (m)
    %    - dhdt_annual: ice thickness change over at each year (m/yr)
    %    - dhdt_monthly: ice thickness change at monthly intervals (m/month)
    %    - years_thickness: timearray on which 'h_annual' is defined (yr)
    %    - lat_sphere: latitude coordinates on sphere
    %    - long_sphere: longitude coordinates on sphere
    %    - X: original X grid coordinates (m)
    %    - Y: original Y grid coordinates (m)
    
    % Note: rhoo and rhoi should be defined as constants
    rhoo = 1000.0; % density of water (kg/m^3)
    rhoi = 917.0;  % density of ice (kg/m^3)
    
    format long
    if strcmp(data_name, 'measureItsLive-GEMB') || strcmp(data_name, 'measureItsLive-GSFC')
        disp("Using ice elevation data from MEaSUREs ITS_LIVE (Nilsson and Gardner, 2024)")
        filename = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/nilsson_et_al/Synthesis_GrIS_Perf_r1920m_1992_2023_ESSD_Compliant_v2.nc';
    
        % Always apply firn correction for measureItsLive data - choose the appropriate model
        if strcmp(data_name, 'measureItsLive-GEMB')
            disp("Applying the GEMB firn model to correct the ice thickness change")
            firnmodel = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/nilsson_et_al/GEMB_FAC_SMB_ALT_GRID.h5'; % FAC model Glacier Energy and Mass Balance FDM (GEMB) version 1.2 (Gardner et al., 2023) 
        elseif strcmp(data_name, 'measureItsLive-GSFC')
            disp("Applying the GSFC firn model to correct the ice thickness change")
            firnmodel = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/nilsson_et_al/GSFC_GrIS_FAC_SMB_ALT_GRID.h5';% FAC model Goddard Space Flight Center FDM version 1.2.1 (Medley et al., 2022).
        else
            error('preprocess_ice_altimetry: Unsupported measureItsLive dataset. Expected measureItsLive-GEMB or measureItsLive-GSFC');
        end
    

        % Read in the variables for the ice elevation data
        %X = ncread(filename, 'x');       % East coordinate (m)
        X = linspace(-644167, 857273, 783); % manually set the X coordinates for now because the original data is not regular (bug)
        Y = ncread(filename, 'y');       % North coordinate (m)
        Time = ncread(filename, 'time'); % Seconds since 1950-01-01
        dh = ncread(filename, 'dh');   % elevation change in meters relative to 2014-01-01, has size of [len(Time) len(Y) len(X)]
        mask_cont = ncread(filename, 'mask_cont'); % ice sheet mask (Zwally et al., 2012)
        rms_sys = ncread(filename, 'rms_sys');
        rms_ran = ncread(filename, 'rms_ran');
    
        % Below is needed (for some reason) to read in the mask_perf variables
        ncid = netcdf.open(filename, 'NOWRITE'); % Open the file
        varid = netcdf.inqVarID(ncid, 'mask_perf_cl01'); % Get the variable ID for mask_perf_cl01
        mask_perf_cl01 = netcdf.getVar(ncid, varid); % Read the variable
        varid = netcdf.inqVarID(ncid, 'mask_perf_cl02'); % Get the variable ID for mask_perf_cl01
        mask_perf_cl02 = netcdf.getVar(ncid, varid); % Read the variable
        netcdf.close(ncid); % Close the file
        
        % Debug statement
        % disp(['Size of mask_cont: ', num2str(size(mask_cont))]);
        % disp(['Size of mask_perf_cl01: ', num2str(size(mask_perf_cl01))]);
        % disp(['Size of mask_perf_cl02: ', num2str(size(mask_perf_cl02))]);
        % disp(['Size of rms_sys: ', num2str(size(rms_sys))]);
        % disp(['Size of rms_ran: ', num2str(size(rms_ran))]);
    
        % Set fill values to NaN in 'dh'
        dh(dh == -32767) = NaN;

        % Apply firn correction
        dfac = ncread(firnmodel, 'dfac');

        % Apply the firn correction to the ice thickness for monthly interval
        dh_corr = dh - dfac;

        % Isolate month-to-month changes
        dhdt_monthly = diff(dh_corr,1,3); % get first difference along the 'time' axis
        dhdt_monthly = dhdt_monthly(:, :, 12:end); % get data from 1993-2023 (i.e. skip the first year)
        Time = Time(12:end);

        % Convert time from seconds since 1950-01-01 to decimal years
        reference_date = datetime(1950, 1, 1);
        Time_datetime = reference_date + seconds(Time);
        Time_years = year(Time_datetime) + (day(Time_datetime, 'dayofyear') - 1) / 365.25;
        
        % Get time vector for the data we'll extract from actual Time data
        years = unique(floor(Time_years)); % get unique years
        years = years(2:end); % remove the first year (incomplete)

    elseif strcmp(data_name, 'DTU2016')
        disp("Using ice elevation data from Khan et al. 2016")
        filename = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/khan_et_al/Greenland_dhdt_mass_1kmgrid.nc';
        
        % Read the variables
        X = ncread(filename, 'X');       % East coordinate (km)
        Y = ncread(filename, 'Y');       % North coordinate (km)
        Time = ncread(filename, 'Time'); % Time in years
        dhdt_water = ncread(filename, 'dhdt'); % Mean monthly elevation change rate in water equivalent height corrected for GIA, elastic uplift and firn compaction in water equivalent height
        
        % Convert dhdt from water equivalent to ice thickness change
        dhdt_monthly = dhdt_water .* (rhoo / rhoi);
    
        % Filter data to end of 2022 (exclude incomplete 2023 data)
        time_mask = Time <= 2023.0; % Include only data up to end of 2022
        Time = Time(time_mask);
        dhdt_monthly = dhdt_monthly(:, :, time_mask);

        % Define years vector for output from actual Time data
        years = unique(floor(Time));
    
        % convert units to meters
        X = X.*1000.0;
        Y = Y.*1000.0;
    
    elseif strcmp(data_name, 'DTU2025')
        disp("Using ice elevation data from DTU2025")
        filename = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/altimetry/khan_et_al_2025/Greenland_dhdt_icevol_1kmgrid_DB.nc';
        fname_firn = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/altimetry/khan_et_al_2025/Greenland_dhdt_firn_1kmgrid_DB.nc';
        
        % Read the variables
        X = ncread(filename, 'x');       % East coordinate (km)
        Y = ncread(filename, 'y');       % North coordinate (km)
        Time = ncread(filename, 'time'); % Time in days since 2003-01-01
        dhdt_monthly = ncread(filename, 'dhdt_vol'); % Mean monthly elevation change rate in ice equivalent height

        %dhdt_firn = ncread(fname_firn, 'dhdt_firn'); % Mean monthly elevation change rate in ice equivalent height

        % Convert time from days since 2003-01-01 to decimal years
        reference_date = datetime(2003, 1, 1);
        Time_datetime = reference_date + days(Time);
        Time = year(Time_datetime) + (day(Time_datetime, 'dayofyear') - 1) / 365.25;
        
        % Filter data to end of 2022 (exclude incomplete 2023 data)
        time_mask = Time <= 2023.0; % Include only data up to end of 2022
        Time = Time(time_mask);
        dhdt_monthly = dhdt_monthly(:, :, time_mask);
        %dhdt_firn = dhdt_firn(:, :, time_mask);
        
        % Get uncorrected dhdt_monthly - comment out for now
        %dhdt_monthly = dhdt_water + dhdt_firn;

        % Define years vector for output from actual Time data
        years = unique(floor(Time));

    elseif strcmp(data_name, 'Buffalo2025-GEMB') || strcmp(data_name, 'Buffalo2025-GSFC') || strcmp(data_name, 'Buffalo2025-IMAU') % Note their data does not resolve monthly data
        disp("Using ice elevation data from Gao et al. 2025")
        filename_GEMB = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/altimetry/Gao_et_al_2025/Annual_rates_grids_EPSG3413/MB_Greenland_GEMB_1994_2020.nc'; %These are the same file but ncks-ed to remove the last year (2021)
        filename_GSFC = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/altimetry/Gao_et_al_2025/Annual_rates_grids_EPSG3413/MB_Greenland_GSFC_1994_2020.nc'; %These are the same file but ncks-ed to remove the last year (2021)
        filename_IMAU = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/altimetry/Gao_et_al_2025/Annual_rates_grids_EPSG3413/MB_Greenland_IMAU_1994_2020.nc';

        % Choose the correct filename based on data_name
        if strcmp(data_name, 'Buffalo2025-GEMB')
            filename = filename_GEMB;
        elseif strcmp(data_name, 'Buffalo2025-GSFC')
            filename = filename_GSFC;
        elseif strcmp(data_name, 'Buffalo2025-IMAU')
            filename = filename_IMAU;
        end

        % Read the variables
        X = ncread(filename, 'x');       % East coordinate (km)
        Y = ncread(filename, 'y');       % North coordinate (km)
        Time = ncread(filename, 'time'); % Time in years
        dhdt_annual = ncread(filename, 'dh'); % Annual elevation change between Sep 1 of 1994 and 2020
        % Define years vector for output from actual Time data
        years = Time; % Time is already in years for this dataset
        
    end
    
    % Display size of the variables to verify (debug statement)
    disp(['Size of X: ', num2str(size(X))]);
    disp(['Size of Y: ', num2str(size(Y))]);
    disp(['Size of Time: ', num2str(size(Time))]);
    if strcmp(data_name, 'measureItsLive') || strcmp(data_name, 'DTU2016')
        disp(['Size of dhdt_monthly: ', num2str(size(dhdt_monthly))]);
    end
    
    % Handle different data structures
    if strcmp(data_name, 'DTU2025')
        % DTU2025 data is already in monthly intervals, just need to group by years
        % Create annual data by averaging monthly data for each year
        unique_years = unique(floor(Time));
        num_years = length(unique_years);
        dhdt_annual = zeros(size(dhdt_monthly, 1), size(dhdt_monthly, 2), num_years);
        
        for i = 1:num_years
            year_mask = floor(Time) == unique_years(i);
            if any(year_mask)
                dhdt_annual(:,:,i) = mean(dhdt_monthly(:,:,year_mask), 3, 'omitnan')*12;
            end
        end
        years = unique_years;
    elseif strcmp(data_name, 'measureItsLive-GEMB') || strcmp(data_name, 'measureItsLive-GSFC') || strcmp(data_name, 'DTU2016')
        % MEaSUREs and DTU2016 data - use original reshaping logic
        num_years = size(dhdt_monthly, 3) / 12;
        dhdt_monthly_reshaped = reshape(dhdt_monthly, size(dhdt_monthly, 1), size(dhdt_monthly, 2), 12, num_years);
        
        % Sum over the 12 months (3rd dimension) to get annual changes ignoring NaNs
        dhdt_annual = sum(dhdt_monthly_reshaped, 3, 'omitnan');
        dhdt_annual = squeeze(dhdt_annual);
    elseif strcmp(data_name, 'Buffalo2025-GEMB') || strcmp(data_name, 'Buffalo2025-GSFC') || strcmp(data_name, 'Buffalo2025-IMAU')
        num_years = size(dhdt_annual, 3);
        % Other than that,do nothing since the data is already in annual intervals
        dhdt_monthly = NaN;
    end
    
    % Calculate the total change across all years, ignoring NaNs
    dhdt_total = sum(dhdt_annual, 3, 'omitnan');

    % Convert years to double precision for polyfit compatibility in get_ice_mass_change function
    years = double(years);

    % Create years_thickness for thickness data (has nt+1 elements)
    % years represents years for dhdt_annual, years_thickness represents years for h_annual
    years_thickness = [years; years(end) + 1]; % Add one more year for the final thickness

    % Calculate ice mass change in Gt
    disp("Calculating ice mass change based on annual data")
    [dice_mass_annual] = get_ice_mass_change(rhoi, X, Y, dhdt_annual, years, false);

    %disp("Calculating ice mass change based on monthly data - Debug")
    %[dice_mass_monthly] = get_ice_mass_change(rhoi, X, Y, dhdt_monthly, Time, false); % check based on the monthly data for the DTU data
    %[dice_mass_monthly] = get_ice_mass_change(rhoi, X, Y, dhdt_monthly, Time_years(2:end), false); % check based on the monthly data for MEaSUREs

    % Transform decimal years to yyyy-mm-dd
    % Extract the integer year
    years_part = floor(Time);
    
    % Calculate the fraction of the year
    % fractions_part = Time - years_part;% Determine if it's a leap year
    % is_leap = leapyear(years_part);
    % days_in_year = 365 + is_leap; 
    
    % Convert the fractional part to days
    % fractional_days = fractions_part .* days_in_year;
    
    % Add the fractional days to the start of the year
    % date = datetime(years_part, 1, 1) + days(fractional_days);
    
    % Create projection structure using EPSG code
    proj_info = projcrs(3413);
    projcrs(3413).ProjectionParameters;

    % If X and Y are vectors of different size, create a meshgrid for the coordinates
    if ((isvector(X) && isvector(Y)) && ~isequal(size(X), size(Y)))
        X = double(X(:)); % Ensure X and Y have the same orientation and floating points
        Y = double(Y(:));
        [X_2d, Y_2d] = meshgrid(X, Y);
    end

    % Perform inverse projection to get latitude and longitude on the WGS84
    % ellipsoid
    [lat_ellipsoid, long_ellipsoid] = projinv(proj_info, X_2d, Y_2d); % put x and y in meters

    % TRANSFORM TO SPHERE for GIA modeling (supports full 3D deformation)
    disp('Transforming to spherical coordinates');
    r_earth = 6371000.; % radius of the earth consistent with the value used in ISSM
    [lat_sphere, long_sphere, dhdt_annual] = ellipsoid_to_sphere(lat_ellipsoid, long_ellipsoid, r_earth, dhdt_annual);

    % Convert NaN values to zero (no change in ice thickness)
    fprintf('Converting NaN values to zero...\n');
    dhdt_annual(isnan(dhdt_annual)) = 0;
    fprintf('After NaN conversion: dhdt_annual range = %.1f to %.1f m/yr\n', ...
        min(dhdt_annual(:)), max(dhdt_annual(:)));
    fprintf('=====================================\n');

    % Calculate total change across years for plotting
    dhdt_total = sum(dhdt_annual, 3, 'omitnan');

    % Create a time series of ice thickness
    % h_annual should have nt+1 time steps (thickness at start of each year)
    h_annual = zeros(size(dhdt_annual, 1), size(dhdt_annual, 2), size(dhdt_annual, 3) + 1);

    % Debug: Check dhdt_annual values
    fprintf('=== DEBUG: dhdt_annual statistics ===\n');
    fprintf('Size of dhdt_annual: %s\n', mat2str(size(dhdt_annual)));

    % Check for NaN values
    nan_count = sum(isnan(dhdt_annual(:)));
    total_count = numel(dhdt_annual);
    fprintf('NaN values: %d (%.2f%% of total)\n', nan_count, 100*nan_count/total_count);

    % Calculate statistics for all values
    fprintf('Min dhdt_annual: %.3f m/yr\n', min(min(dhdt_annual(:))));
    fprintf('Max dhdt_annual: %.3f m/yr\n', max(max(dhdt_annual(:))));
    fprintf('Mean dhdt_annual: %.3f m/yr\n', mean(mean(dhdt_annual(:))));
    fprintf('Std dhdt_annual: %.3f m/yr\n', std(std(dhdt_annual(:))));

    % Check cumulative loss
    cumulative_loss = sum(dhdt_annual, 3, 'omitnan');
    fprintf('Min cumulative loss: %.3f m\n', min(cumulative_loss(:), [], 'omitnan'));
    fprintf('Max cumulative loss: %.3f m\n', max(cumulative_loss(:), [], 'omitnan'));
    fprintf('Mean cumulative loss: %.3f m\n', mean(cumulative_loss(:), 'omitnan'));
    fprintf('Std cumulative loss: %.3f m\n', std(cumulative_loss(:), 'omitnan'));

    % Check how many pixels have cumulative loss > 1500m
    pixels_over_1500m = sum(cumulative_loss(:) > 1500);
    total_pixels = numel(cumulative_loss);
    fprintf('Pixels with cumulative loss > 1500m: %d (%.2f%% of total)\n', pixels_over_1500m, 100*pixels_over_1500m/total_pixels);

    % Create time series of ice thickness using forward integration
    % Start with a realistic initial thickness and integrate forward
    fprintf('=== Reconstructing ice thickness using forward integration approach ===\n');

    % Set initial thickness (start of first year)
    initial_thickness = 1500; % Start with 1500m as initial thickness
    h_annual(:,:,1) = initial_thickness;

    % Integrate forward: h(t+1) = h(t) + dhdt(t)
    for t = 1:size(dhdt_annual,3)

        % get indices of dhdt_annual values that are NaN
        nan_indices_dhdt = find(isnan(dhdt_annual(:,:,t)));
        if ~isempty(nan_indices_dhdt)
            fprintf('dhdt_annual is NaN at %d indices\n', length(nan_indices_dhdt));
        end

        nan_indices_h = find(isnan(h_annual(:,:,t)));
        if ~isempty(nan_indices_h)
            fprintf('h_annual is NaN at %d indices\n', length(nan_indices_h));
        end

        % Update thickness
        h_annual(:,:,t+1) = h_annual(:,:,t) + dhdt_annual(:,:,t);

        fprintf('At step %d: dhdt_annual range = %.3f to %.3f m/yr\n', t, min(min(dhdt_annual(:,:,t))), max(max(dhdt_annual(:,:,t))));
        fprintf('At step %d: h_annual range = %.3f to %.3f m\n', t+1, min(min(h_annual(:,:,t+1))), max(max(h_annual(:,:,t+1))));
    end

    h_annual = permute(h_annual, [2, 1, 3]); % flip the x and y axes
    dhdt_annual = permute(dhdt_annual, [2, 1, 3]); % flip the x and y axes
    dhdt_total = permute(dhdt_total, [2, 1]); % flip the x and y axes

    if strcmp(data_name, 'measureItsLive') || strcmp(data_name, 'DTU2016')
        dhdt_monthly = permute(dhdt_monthly, [2, 1, 3]); % flip the x and y axes
    end

    disp('====================================');

    % Optional plotting
    if plot_altimetry
        base_fig_num = 100;
        disp('Plotting figures...')
        num_colors = 100; 
        Bl2white2Rd_cmap = [linspace(0, 1, num_colors)', linspace(0, 1, num_colors)', ones(num_colors, 1);  % Blue to white
            ones(num_colors, 1), linspace(1, 0, num_colors)', linspace(1, 0, num_colors)']; % White to red

        % Plot the ice thickness change from the original data
        data =  dhdt_total; %h_annual(:,:,end) - h_annual(:,:,1);
        %figure(base_fig_num + 1);
        %clf; % Clear the figure
        %imagesc(X, Y, data);
        %set(gca, 'YDir', 'normal'); % Correct orientation if needed
        %xlabel('East (km)','FontSize',14);
        %ylabel('North (km)','FontSize',14);
        %title(sprintf('Total ice elevation change on x-y plane between %d-%d (m) %s', years(1), years(end), data_name), 'FontSize', 14);
        %colorbarHandle = colorbar;  % Create the colorbar and get its handle
        %colorbarHandle.FontSize = 14;
        %colormap(flipud(Bl2white2Rd_cmap));
        %caxis([-50 50]);
    
        figure(base_fig_num + 2)
        p=pcolor(long_sphere,lat_sphere, data);
        shading flat;
        set(gca, 'YDir', 'normal'); % Correct orientation if needed
        set(p, 'AlphaData', ~isnan(data));  % Make NaN regions transparent
        set(gca, 'Color', 'white');  
        xlabel('Longitude (degrees)','FontSize',14);
        ylabel('Latitude (degrees)','FontSize',14);
        title(sprintf('Total ice elevation change between %d-%d (m) %s', years(1), years(end), data_name), 'FontSize', 14);
        colorbarHandle = colorbar;  % Create the colorbar and get its handle
        colorbarHandle.FontSize = 14;
        colormap(flipud(Bl2white2Rd_cmap));
        %caxis([-50 50]);
    
        % plot mean elevation change m/yr
        figure(base_fig_num + 3)
        p=pcolor(long_sphere,lat_sphere, data/num_years);
        shading flat;
        set(gca, 'YDir', 'normal'); % Correct orientation if needed
        set(p, 'AlphaData', ~isnan(data));  % Make NaN regions transparent
        set(gca, 'Color', 'white');  
        xlabel('Longitude (degrees)','FontSize',14);
        ylabel('Latitude (degrees)','FontSize',14);
        title(sprintf('Mean ice elevation change between %d-%d (m/yr) %s', years(1), years(end), data_name), 'FontSize', 14);
        colorbarHandle = colorbar;  % Create the colorbar and get its handle
        colorbarHandle.FontSize = 14;
        colormap(flipud(Bl2white2Rd_cmap));
        caxis([-1 1]);
    
        % Plot an annual ice elevation change
        data = dhdt_annual(:,:,end);
        figure(base_fig_num + 4)
        p=pcolor(long_sphere,lat_sphere, data);
        shading flat;
        % Make NaN regions appear as white by setting AlphaData
        set(p, 'AlphaData', ~isnan(data));  % Make NaN regions transparent
        set(gca, 'Color', 'white');          % Set background color to white
        xlabel('Longitude (degrees)','FontSize',14);
        ylabel('Latitude (degrees)','FontSize',14);
        title(sprintf('Ice elevation change between %d-%d (m/yr) %s', years(end-1), years(end), data_name), 'FontSize', 14);
        colorbarHandle = colorbar;  % Create the colorbar and get its handle
        colorbarHandle.FontSize = 14;
        colormap(flipud(Bl2white2Rd_cmap));
        caxis([-1 1]);
    end
end 