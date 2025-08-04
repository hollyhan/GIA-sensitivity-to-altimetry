function [dfac_annual, dfac_monthly, years, lat_sphere, long_sphere] = preprocess_firn_model(firnmodel_name)
    % preprocess_firn_model.m
    % Holly Han (created: July 25th, 2025; Last edited: July 25th, 2025).
    % Preprocesses firn air compaction data from different models.
    %
    % Inputs:
    %   - firnmodel_name: Model name ('GEMB' or 'GSFC')
    % Outputs:
    %    - dfac_annual: firn air compaction at annual time interval (m)
    %    - dfac_monthly: firn air compaction at monthly intervals (m/month)
    %    - years: timearray on which 'dfac_annual' is defined (yr)
    %    - lat_sphere: latitude coordinates on sphere
    %    - long_sphere: longitude coordinates on sphere
    
    format long
    
    % Define file paths for firn models
    if strcmp(firnmodel_name, 'GEMB')
        disp("Loading GEMB firn model data")
        firnmodel = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/nilsson_et_al/GEMB_FAC_SMB_ALT_GRID.h5'; % FAC model Glacier Energy and Mass Balance FDM (GEMB) version 1.2 (Gardner et al., 2023) 
    elseif strcmp(firnmodel_name, 'GSFC')
        disp("Loading GSFC firn model data")
        firnmodel = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/nilsson_et_al/GSFC_GrIS_FAC_SMB_ALT_GRID.h5';% FAC model Goddard Space Flight Center FDM version 1.2.1 (Medley et al., 2022).
    elseif strcmp(firnmodel_name, 'RACMO2.3p2')
        disp("Loading RACMO2.3p2 firn model data")
        firnmodel = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/altimetry/khan_et_al_2025/Greenland_dhdt_firn_1kmgrid_DB.nc'; % FAC model RACMO2.3p2 (Khan et al., 2025)
    end
    
    % Read in variables for firn correction
    if strcmp(firnmodel_name, 'GEMB')
        % GEMB model variables (HDF5 format)
        X = ncread(firnmodel, 'X');       % East coordinate (m)
        Y = ncread(firnmodel, 'Y');       % North coordinate (m)
        Time = ncread(firnmodel, 'time'); % Time in decimal years
        dfac = ncread(firnmodel, 'dfac'); % Firn air content (m)
        years = 1993:2023;
    elseif strcmp(firnmodel_name, 'GSFC')
        % GSFC model variables (HDF5 format)
        X = ncread(firnmodel, 'X');       % East coordinate (m)
        Y = ncread(firnmodel, 'Y');       % North coordinate (m)
        Time = ncread(firnmodel, 'time'); % Time in decimal years
        dfac = ncread(firnmodel, 'dfac');  % Firn air content (m)
        years = 1993:2023;
    elseif strcmp(firnmodel_name, 'RACMO2.3p2')
        % DTU2025 model variables (NetCDF format)
        X = ncread(firnmodel, 'x');       % East coordinate (m)
        Y = ncread(firnmodel, 'y');       % North coordinate (m)
        Time = ncread(firnmodel, 'time'); % Time in decimal years
        dfac = ncread(firnmodel, 'dhdt_firn'); % Firn air content (m)
        
        % Convert time from days since 2003-01-01 to decimal years
        reference_date = datetime(2003, 1, 1);
        Time_datetime = reference_date + days(Time);
        Time = year(Time_datetime) + (day(Time_datetime, 'dayofyear') - 1) / 365.25;
        
        % Filter data to end of 2022 (exclude incomplete 2023 data)
        time_mask = Time <= 2023.0; % Include only data up to end of 2022
        Time = Time(time_mask);
        dfac = dfac(:, :, time_mask);
        years = 2003:2022;
    end
    
    % Store monthly data for output
    if strcmp(firnmodel_name, 'RACMO2.3p2')
        % RACMO data is already mean monthly values, no differencing needed
        dfac_monthly = dfac;
    else
        % GEMB and GSFC data are cumulative, need differencing
        dfac_monthly = diff(dfac, 1, 3);
        dfac_monthly = dfac_monthly(:,:,12:end); %get data from 1993 to 2023
    end
    
    % Calculate annual firn air compaction by averaging monthly data
    num_years = size(dfac_monthly, 3) / 12;
    dfac_monthly_reshaped = reshape(dfac_monthly, size(dfac_monthly, 1), size(dfac_monthly, 2), 12, num_years);
    
    % Sum over the 12 months (3rd dimension) to get annual changes ignoring NaNs
    dfac_annual = sum(dfac_monthly_reshaped, 3, 'omitnan');
    dfac_annual = squeeze(dfac_annual);
    
    % Create projection structure using EPSG code
    proj_info = projcrs(3413);
    
    % If X and Y are vectors of different size, create a meshgrid for the coordinates
    if ((isvector(X) && isvector(Y)) && ~isequal(size(X), size(Y)))
        X = double(X(:)); % Ensure X and Y have the same orientation and floating points
        Y = double(Y(:));
        [X_2d, Y_2d] = meshgrid(X, Y);
    else
        X_2d = X;
        Y_2d = Y;
    end
    
    % Perform inverse projection to get latitude and longitude on the WGS84 ellipsoid
    [lat_ellipsoid, long_ellipsoid] = projinv(proj_info, X_2d, Y_2d);
    
    % TRANSFORM TO SPHERE for GIA modeling (supports full 3D deformation)
    disp('Transforming to spherical coordinates for GIA modeling (supports 3D deformation)');
    r_earth = 6371000; % Earth radius in meters
    [lat_sphere, long_sphere, dfac_annual] = ellipsoid_to_sphere(lat_ellipsoid, long_ellipsoid, r_earth, dfac_annual, false);
    
    % Optional plotting
    if true % Set to true to enable plotting
        base_fig_num = 200;
        disp('Plotting firn model figures...')
        num_colors = 100; 
        Bl2white2Rd_cmap = [linspace(0, 1, num_colors)', linspace(0, 1, num_colors)', ones(num_colors, 1);  % Blue to white
            ones(num_colors, 1), linspace(1, 0, num_colors)', linspace(1, 0, num_colors)']; % White to red

        % Plot the total firn air compaction change
        data = sum(dfac_annual, 3, 'omitnan');
        figure(base_fig_num + 1);
        clf; % Clear the figure
        % Debug: check dimensions
        disp(['Size of X: ', num2str(size(X))]);
        disp(['Size of Y: ', num2str(size(Y))]);
        disp(['Size of data: ', num2str(size(data))]);
        
        % Handle different coordinate formats
        if isvector(X) && isvector(Y)
            % X and Y are vectors, create meshgrid for plotting
            [X_plot, Y_plot] = meshgrid(X, Y);
            pcolor(X_plot, Y_plot, data');  % Note: transpose data for meshgrid convention
        else
            % X and Y are already 2D arrays
            pcolor(X, Y, data);
        end
        shading flat;
        set(gca, 'YDir', 'normal'); % Correct orientation if needed
        xlabel('East (m)','FontSize',14);
        ylabel('North (m)','FontSize',14);
        title(sprintf('Total firn air compaction change %s model', firnmodel_name), 'FontSize', 14);
        colorbarHandle = colorbar;  % Create the colorbar and get its handle
        colorbarHandle.FontSize = 14;
        colormap(flipud(Bl2white2Rd_cmap));
        %caxis([-10 10]);
    
        figure(base_fig_num + 2)
        % Debug: check dimensions for spherical plotting
        disp(['Size of long_sphere: ', num2str(size(long_sphere))]);
        disp(['Size of lat_sphere: ', num2str(size(lat_sphere))]);
        disp(['Size of data for spherical: ', num2str(size(data))]);
        
        % Handle different coordinate formats for spherical plotting
        if isvector(long_sphere) && isvector(lat_sphere)
            % Coordinates are vectors, create meshgrid for plotting
            [long_plot, lat_plot] = meshgrid(long_sphere, lat_sphere);
            p=pcolor(long_plot, lat_plot, data');  % Note: transpose data for meshgrid convention
        else
            % Coordinates are already 2D arrays - check if dimensions match
            if ~isequal(size(long_sphere), size(data))
                % Transpose data to match coordinate dimensions
                p=pcolor(long_sphere, lat_sphere, data');
            else
                p=pcolor(long_sphere, lat_sphere, data);
            end
        end
        shading flat;
        set(gca, 'YDir', 'normal'); % Correct orientation if needed
        set(p, 'AlphaData', ~isnan(data));  % Make NaN regions transparent
        set(gca, 'Color', 'white');  
        xlabel('Longitude (degrees)','FontSize',14);
        ylabel('Latitude (degrees)','FontSize',14);
        title(sprintf('Total firn air compaction change %s model (spherical)', firnmodel_name), 'FontSize', 14);
        colorbarHandle = colorbar;  % Create the colorbar and get its handle
        colorbarHandle.FontSize = 14;
        colormap(flipud(Bl2white2Rd_cmap));
        %caxis([-10 10]);
    
        % Plot mean annual firn air compaction
        figure(base_fig_num + 3)
        % Handle different coordinate formats for spherical plotting
        if isvector(long_sphere) && isvector(lat_sphere)
            % Coordinates are vectors, create meshgrid for plotting
            [long_plot, lat_plot] = meshgrid(long_sphere, lat_sphere);
            p=pcolor(long_plot, lat_plot, (data/size(dfac_annual, 3))');  % Note: transpose data for meshgrid convention
        else
            % Coordinates are already 2D arrays - check if dimensions match
            if ~isequal(size(long_sphere), size(data))
                % Transpose data to match coordinate dimensions
                p=pcolor(long_sphere, lat_sphere, (data/size(dfac_annual, 3))');
            else
                p=pcolor(long_sphere, lat_sphere, data/size(dfac_annual, 3));
            end
        end
        shading flat;
        set(gca, 'YDir', 'normal'); % Correct orientation if needed
        set(p, 'AlphaData', ~isnan(data));  % Make NaN regions transparent
        set(gca, 'Color', 'white');  
        xlabel('Longitude (degrees)','FontSize',14);
        ylabel('Latitude (degrees)','FontSize',14);
        title(sprintf('Mean annual firn air compaction %s model', firnmodel_name), 'FontSize', 14);
        colorbarHandle = colorbar;  % Create the colorbar and get its handle
        colorbarHandle.FontSize = 14;
        colormap(flipud(Bl2white2Rd_cmap));
        %caxis([-0.5 0.5]);
    
        % Plot an annual firn air compaction
        data_annual = dfac_annual(:,:,end);
        figure(base_fig_num + 4)
        % Handle different coordinate formats for spherical plotting
        if isvector(long_sphere) && isvector(lat_sphere)
            % Coordinates are vectors, create meshgrid for plotting
            [long_plot, lat_plot] = meshgrid(long_sphere, lat_sphere);
            p=pcolor(long_plot, lat_plot, data_annual');  % Note: transpose data for meshgrid convention
        else
            % Coordinates are already 2D arrays - check if dimensions match
            if ~isequal(size(long_sphere), size(data_annual))
                % Transpose data to match coordinate dimensions
                p=pcolor(long_sphere, lat_sphere, data_annual');
            else
                p=pcolor(long_sphere, lat_sphere, data_annual);
            end
        end
        shading flat;
        % Make NaN regions appear as white by setting AlphaData
        set(p, 'AlphaData', ~isnan(data_annual));  % Make NaN regions transparent
        set(gca, 'Color', 'white');          % Set background color to white
        xlabel('Longitude (degrees)','FontSize',14);
        ylabel('Latitude (degrees)','FontSize',14);
        title(sprintf('Annual firn air compaction %s model in year %d', firnmodel_name, years(end)), 'FontSize', 14);
        colorbarHandle = colorbar;  % Create the colorbar and get its handle
        colorbarHandle.FontSize = 14;
        colormap(flipud(Bl2white2Rd_cmap));
        %caxis([-0.5 0.5]);
    end
end 