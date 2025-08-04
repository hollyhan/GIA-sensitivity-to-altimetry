function [h_annual, dhdt_annual, dhdt_monthly, years, lat_sphere, long_sphere] = preprocess_ice_altimetry(data_name)
    % preprocess_ice_altimetry.m
    % Holly Han (created: July 25th, 2025; Last edited: July 25th, 2025).
    % Preprocesses ice thickness and elevation change data from altimetry.
    %
    % Inputs:
    %   - data_name: Dataset name ('measureItsLive', 'DTU2016', or 'DTU2024')
    % Outputs:
    %    - h_annual: ice thickness at annual time interval (m)
    %    - dhdt_annual: ice thickness change over at each year (m/yr)
    %    - dhdt_monthly: ice thickness change at monthly intervals (m/month)
    %    - years: timearray on which 'h_annual' and 'dhdt_annual' are defined (yr)
    %    - lat_sphere: latitude coordinates on sphere
    %    - long_sphere: longitude coordinates on sphere
    
    % Note: rhoo and rhoi should be defined as constants
    rhoo = 1000.0; % density of water (kg/m^3)
    rhoi = 917.0;  % density of ice (kg/m^3)
    
    format long
    if strcmp(data_name, 'measureItsLive')
        disp("Using ice elevation data from MEaSUREs ITS_LIVE (Nilsson and Gardner, 2024)")
        filename = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/nilsson_et_al/Synthesis_GrIS_Perf_r1920m_1992_2023_ESSD_Compliant_v2.nc';
    
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
        
        % Isolate month-to-month changes
        dhdt_monthly = diff(dh,1,3); % get first difference along the 'time' axis
        dhdt_monthly = dhdt_monthly(:, :, 12:end); % get data from 1993-2023
    
        % Get time vector for the data we'll extract
        years = 1993:2023;
        
        % Note: Firn correction is now handled separately by preprocess_firn_model.m
        % To apply firn correction, call preprocess_firn_model() and subtract dfac_annual
        % from dhdt_annual after processing
    
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
    
        % Define years vector for output
        years = 2003:2022;
        dhdt_monthly = dhdt_monthly(:, :, 1:240);
    
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
        dhdt_water = ncread(filename, 'dhdt_vol'); % Mean monthly elevation change rate in ice equivalent height
        dhdt_firn = ncread(fname_firn, 'dhdt_firn'); % Mean monthly elevation change rate in ice equivalent height

        % Convert time from days since 2003-01-01 to decimal years
        reference_date = datetime(2003, 1, 1);
        Time_datetime = reference_date + days(Time);
        Time = year(Time_datetime) + (day(Time_datetime, 'dayofyear') - 1) / 365.25;
        
        % Filter data to end of 2022 (exclude incomplete 2023 data)
        time_mask = Time <= 2023.0; % Include only data up to end of 2022
        Time = Time(time_mask);
        dhdt_water = dhdt_water(:, :, time_mask);
        dhdt_firn = dhdt_firn(:, :, time_mask);
        
        % Get uncorrected dhdt_monthly
        dhdt_monthly = dhdt_water + dhdt_firn;

        % Define years vector for output (2003-2022)
        years = 2003:2022;
    end
    
    % Display size of the variables to verify (debug statement)
    disp(['Size of X: ', num2str(size(X))]);
    disp(['Size of Y: ', num2str(size(Y))]);
    disp(['Size of Time: ', num2str(size(Time))]);
    disp(['Size of dhdt_monthly: ', num2str(size(dhdt_monthly))]);
    
    % Handle different data structures
    if strcmp(data_name, 'DTU2025')
        % DTU2025 data is already monthly, just need to group by years
        % Create annual data by averaging monthly data for each year
        unique_years = unique(floor(Time));
        num_years = length(unique_years);
        dhdt_annual = zeros(size(dhdt_monthly, 1), size(dhdt_monthly, 2), num_years);
        
        for i = 1:num_years
            year_mask = floor(Time) == unique_years(i);
            if any(year_mask)
                dhdt_annual(:,:,i) = mean(dhdt_monthly(:,:,year_mask), 3, 'omitnan');
            end
        end
        years = unique_years;
    else
        % MEaSUREs and DTU2016 data - use original reshaping logic
        num_years = size(dhdt_monthly, 3) / 12;
        dhdt_monthly_reshaped = reshape(dhdt_monthly, size(dhdt_monthly, 1), size(dhdt_monthly, 2), 12, num_years);
        
        % Sum over the 12 months (3rd dimension) to get annual changes ignoring NaNs
        dhdt_annual = sum(dhdt_monthly_reshaped, 3, 'omitnan');
        dhdt_annual = squeeze(dhdt_annual);
    end
    
    % Calculate the total change across all years, ignoring NaNs
    dhdt_total = sum(dhdt_annual, 3, 'omitnan');
    
    % Note: Ice mass change calculations removed - implement get_ice_mass_change function if needed
    
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
    disp('Transforming to spherical coordinates for GIA modeling (supports 3D deformation)');
    r_earth = 6371000.; % radius of the earth consistent with the value used in ISSM
    [lat_sphere, long_sphere, dhdt_annual] = ellipsoid_to_sphere(lat_ellipsoid, long_ellipsoid, r_earth, dhdt_annual, true);

    % Calculate total change across years for plotting
    dhdt_total = sum(dhdt_annual, 3, 'omitnan');

    % Create a time series of ice thickness
    h_annual = zeros(size(dhdt_annual));
    h_annual(:,:,end) = 1500; % mean average Greenland ice thickness
    for t = length(years)-1 : -1 : 1
        h_annual(:,:,t) = h_annual(:,:,t+1) - dhdt_annual(:,:,t); % Add dh to move backwards in time
    end

    % Optional plotting (set to true to enable)
    plot_altimetry = true;
    if plot_altimetry
        base_fig_num = 100;
        disp('Plotting figures...')
        num_colors = 100; 
        Bl2white2Rd_cmap = [linspace(0, 1, num_colors)', linspace(0, 1, num_colors)', ones(num_colors, 1);  % Blue to white
            ones(num_colors, 1), linspace(1, 0, num_colors)', linspace(1, 0, num_colors)']; % White to red

        % Plot the ice thickness change from the original data
        data =  dhdt_total'; %h_annual(:,:,end)' - h_annual(:,:,1)';
        figure(base_fig_num + 1);
        clf; % Clear the figure
        imagesc(X, Y, data);
        set(gca, 'YDir', 'normal'); % Correct orientation if needed
        xlabel('East (km)','FontSize',14);
        ylabel('North (km)','FontSize',14);
        title(sprintf('Total ice elevation change on x-y plane between %d-%d', years(1), years(end)), 'FontSize', 14);
        colorbarHandle = colorbar;  % Create the colorbar and get its handle
        colorbarHandle.FontSize = 14;
        colormap(flipud(Bl2white2Rd_cmap));
        %caxis([-50 50]);
    
        figure(base_fig_num + 2)
        p=pcolor(long_sphere,lat_sphere, data);
        shading flat;
        set(gca, 'YDir', 'normal'); % Correct orientation if needed
        set(p, 'AlphaData', ~isnan(data));  % Make NaN regions transparent
        set(gca, 'Color', 'white');  
        xlabel('Longitude (degrees)','FontSize',14);
        ylabel('Latitude (degrees)','FontSize',14);
        title(sprintf('Total ice elevation change (spherical) between %d-%d', years(1), years(end)), 'FontSize', 14);
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
        title(sprintf('Mean ice elevation change (spherical) between %d-%d', years(1), years(end)), 'FontSize', 14);
        colorbarHandle = colorbar;  % Create the colorbar and get its handle
        colorbarHandle.FontSize = 14;
        colormap(flipud(Bl2white2Rd_cmap));
        caxis([-1 1]);
    
        % Plot an annual ice elevation change
        data = dhdt_annual(:,:,end)';
        figure(base_fig_num + 4)
        p=pcolor(long_sphere,lat_sphere, data);
        shading flat;
        % Make NaN regions appear as white by setting AlphaData
        set(p, 'AlphaData', ~isnan(data));  % Make NaN regions transparent
        set(gca, 'Color', 'white');          % Set background color to white
        xlabel('Longitude (degrees)','FontSize',14);
        ylabel('Latitude (degrees)','FontSize',14);
        title(sprintf('Annual ice elevation change (spherical) between %d-%d', years(end-1), years(end)), 'FontSize', 14);
        colorbarHandle = colorbar;  % Create the colorbar and get its handle
        colorbarHandle.FontSize = 14;
        colormap(flipud(Bl2white2Rd_cmap));
        caxis([-1 1]);
    end
end 