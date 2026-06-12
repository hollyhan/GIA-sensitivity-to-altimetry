function [h_annual, dhdt_annual, dhdt_monthly, years_thickness, lat_ellipsoid, long_ellipsoid, X, Y, X_2d, Y_2d] = preprocess_ice_altimetry(data_name, plot_altimetry)
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
    %    - lat_ellipsoid: latitude coordinates on ellipsoid
    %    - long_ellipsoid: longitude coordinates on ellipsoid
    %    - X: original X grid coordinates (m)
    %    - Y: original Y grid coordinates (m)
    %    - X_2d: original X grid coordinates in meshgrid format (m)
    %    - Y_2d: original Y grid coordinates in meshgrid format (m)
    
    % Note: rhoo and rhoi should be defined as constants
    rhoo = 1000.0; % density of water (kg/m^3)
    rhoi = 917.0;  % density of ice (kg/m^3)
    
    format long
    if strcmp(data_name, 'measureItsLive-GEMB') || strcmp(data_name, 'measureItsLive-GSFC')
        disp("Using ice elevation data from MEaSUREs ITS_LIVE (Nilsson and Gardner, 2024)")
        filename = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/nilsson_et_al/Greenland_G1920V01_IceSheetGlacierIceHeight.nc';

        % Always apply firn correction for measureItsLive data - choose the appropriate model
        if strcmp(data_name, 'measureItsLive-GEMB')
            disp("Applying the GEMB firn model to correct the ice thickness change")
            dfac_var = 'dfac_gemb';
        elseif strcmp(data_name, 'measureItsLive-GSFC')
            disp("Applying the GSFC firn model to correct the ice thickness change")
            dfac_var = 'dfac_gsfc';
        else
            error('preprocess_ice_altimetry: Unsupported MEaSUREs dataset.');
        end

        % Read in the variables for the ice elevation data
        X = ncread(filename, 'x');       % East coordinate (m)
        Y = ncread(filename, 'y');       % North coordinate (m)
        Time = ncread(filename, 'time'); % days since 1992-01-15 00:00:00

        dh = ncread(filename, 'dh');   % elevation change in meters relative to DEM, w.r.t.  2014-01-01, has size of [len(Time) len(Y) len(X)]
        dfac = ncread(filename, dfac_var);
        mask = ncread(filename, 'mask'); % ice sheet mask (Zwally et al., 2012)
        rms = ncread(filename, 'rms');
        
        % Mask out non-ice regions
        ice_mask = mask == 1 ;%| mask == 2;
        dh(~ice_mask) = NaN;
        dfac(~ice_mask) = NaN;

        % Set fill values to NaN
        dh(dh == -32767) = NaN;
        dfac(dfac == -32767) = NaN;
        rms(rms == -32767) = NaN;

        % Apply the firn correction to the ice thickness for monthly interval
        dh_corr = dh - dfac;

        %% Regression-based mass-rate estimate for MEaSUREs
        % This is for comparison with the publication-style mass budget.
        % It fits dh_corr(t) = a*t + b at each pixel and converts slope a
        % from m/yr to Gt/yr.

        fprintf('\n=== DEBUG MEaSUREs raw fields ===\n');
        fprintf('dh range: %.3f to %.3f m\n', min(dh(:),[],'omitnan'), max(dh(:),[],'omitnan'));
        fprintf('dfac range: %.3f to %.3f m\n', min(dfac(:),[],'omitnan'), max(dfac(:),[],'omitnan'));
        fprintf('dh_corr range: %.3f to %.3f m\n', min(dh_corr(:),[],'omitnan'), max(dh_corr(:),[],'omitnan'));

        ddh = diff(dh,1,3);
        ddfac = diff(dfac,1,3);
        ddhcorr = diff(dh_corr,1,3);

        fprintf('monthly diff dh range: %.3f to %.3f m/month\n', min(ddh(:),[],'omitnan'), max(ddh(:),[],'omitnan'));
        fprintf('monthly diff dfac range: %.3f to %.3f m/month\n', min(ddfac(:),[],'omitnan'), max(ddfac(:),[],'omitnan'));
        fprintf('monthly diff dh_corr range: %.3f to %.3f m/month\n', min(ddhcorr(:),[],'omitnan'), max(ddhcorr(:),[],'omitnan'));
        fprintf('=================================\n\n');

        % Isolate month-to-month changes
        dhdt_monthly = diff(dh_corr, 1, 3); % 383 monthly intervals
        Time = Time(2:end);                 % align Time with diff output

        % Quality control: remove unphysical monthly jumps
        jump_thresh = 10; % m/month
        bad_jump = abs(dhdt_monthly) > jump_thresh;

        fprintf('MEaSUREs QC: monthly jumps > %.1f m/month: %d / %d = %.4f%%\n', ...
            jump_thresh, sum(bad_jump(:)), numel(bad_jump), ...
            100 * sum(bad_jump(:)) / numel(bad_jump));

        dhdt_monthly(bad_jump) = NaN;

        % Trim to start from the second year
        dhdt_monthly = dhdt_monthly(:, :, 12:end);
        Time = Time(12:end);

        % Convert time from days since 1992-01-15 to decimal years
        reference_date = datetime(1992, 1, 15, 0, 0, 0);
        Time_datetime = reference_date + days(Time);
        Time_years = year(Time_datetime) + (day(Time_datetime, 'dayofyear') - 1) / 365.25;
        
        % Get time vector for the data we'll extract from actual Time data
        years = unique(floor(Time_years)); % get unique years

        debug_mask = false;
        if debug_mask
            disp('---Debug statement start----');
            dh_raw = ncread(filename, 'dh');
            dfac_raw = ncread(filename, dfac_var);
            nt = size(dh_raw, 3);

            mask3_gris = repmat(mask == 1, [1 1 nt]);
            mask3_pg   = repmat(mask == 2, [1 1 nt]);
            mask3_all  = repmat(mask == 1 | mask == 2, [1 1 nt]);

            dh_gris = dh_raw; dfac_gris = dfac_raw;
            dh_pg   = dh_raw; dfac_pg   = dfac_raw;
            dh_all  = dh_raw; dfac_all  = dfac_raw;

            dh_gris(~mask3_gris) = NaN;
            dfac_gris(~mask3_gris) = NaN;

            dh_pg(~mask3_pg) = NaN;
            dfac_pg(~mask3_pg) = NaN;

            dh_all(~mask3_all) = NaN;
            dfac_all(~mask3_all) = NaN;

            dhdt_gris = diff(dh_gris - dfac_gris, 1, 3);
            dhdt_pg   = diff(dh_pg   - dfac_pg,   1, 3);
            dhdt_all  = diff(dh_all  - dfac_all,  1, 3);

            % Total change over full differenced period
            dhtot_gris = sum(dhdt_gris, 3, 'omitnan');
            dhtot_pg   = sum(dhdt_pg,   3, 'omitnan');
            dhtot_all  = sum(dhdt_all,  3, 'omitnan');

            % Plot
            figure('Color','w');

            subplot(1,3,1)
            imagesc(X, Y, dhtot_gris');
            set(gca,'YDir','normal')
            axis equal tight
            colorbar
            caxis([-50 50])
            title('GrIS only')
            xlabel('X (m)')
            ylabel('Y (m)')

            subplot(1,3,2)
            imagesc(X, Y, dhtot_pg');
            set(gca,'YDir','normal')
            axis equal tight
            colorbar
            caxis([-50 50])
            title('Peripheral glaciers only')
            xlabel('X (m)')
            ylabel('Y (m)')

            subplot(1,3,3)
            imagesc(X, Y, dhtot_all');
            set(gca,'YDir','normal')
            axis equal tight
            colorbar
            caxis([-50 50])
            title('GrIS + PG')
            xlabel('X (m)')
            ylabel('Y (m)')

            sgtitle(sprintf('%s total corrected elevation change', data_name))

            [dice_mass_gris, dice_mass_cum_gris] = get_ice_mass_change(rhoi, X, Y, dhdt_gris, Time(2:end), false);
            [dice_mass_pg, dice_mass_cum_pg] = get_ice_mass_change(rhoi, X, Y, dhdt_pg, Time(2:end), false);
            [dice_mass_all, dice_mass_cum_all] = get_ice_mass_change(rhoi, X, Y, dhdt_all, Time(2:end), false);

            fprintf('\n==============================\n');
            fprintf('GrIS only total: %.1f Gt\n', dice_mass_cum_gris(end));
            fprintf('PG only total: %.1f Gt\n', dice_mass_cum_pg(end));
            fprintf('GrIS + PG total: %.1f Gt\n', dice_mass_cum_all(end));
            fprintf('Check sum (GrIS + PG): %.1f Gt\n', ...
                dice_mass_cum_gris(end) + dice_mass_cum_pg(end));
            fprintf('==============================\n');
            disp('---Debug statement end----');
        end

        %% Regression-based mass-rate estimate for MEaSUREs
        % This is for comparison with the publication-style mass budget.
        % It fits dh_corr(t) = a*t + b at each pixel and converts slope a
        % from m/yr to Gt/yr.

        run_regression_mass_check = false;
        if run_regression_mass_check
            Time_raw_measures = ncread(filename, 'time');

            reference_date_reg = datetime(1992, 1, 15, 0, 0, 0);
            Time_datetime_reg = reference_date_reg + days(Time_raw_measures);
            Time_years_reg = year(Time_datetime_reg) + ...
                (day(Time_datetime_reg, 'dayofyear') - 1) / 365.25;

            dx_reg = abs(diff(X(1:2)));
            dy_reg = abs(diff(Y(1:2)));
            cell_area_reg = double(dx_reg * dy_reg);

            dhdt_regression = NaN(size(dh_corr,1), size(dh_corr,2));

            fprintf('\nRunning MEaSUREs regression-based mass check...\n');
            for ii = 1:size(dh_corr,1)
                for jj = 1:size(dh_corr,2)
                    ts = squeeze(dh_corr(ii,jj,:));
                    valid = ~isnan(ts);

                    if sum(valid) > 24
                        x = Time_years_reg(valid);
                        y = ts(valid);

                        x = x(:);
                        y = y(:);

                        p = polyfit(x, y, 1);
                        dhdt_regression(ii,jj) = p(1); % m/yr
                    end
                end
            end

            volume_rate_reg = dhdt_regression * cell_area_reg; % m^3/yr
            mass_rate_reg = volume_rate_reg * rhoi * 1.0e-12;  % Gt/yr

            total_mass_rate_reg = sum(mass_rate_reg(:), 'omitnan');

            time_span_reg = Time_years_reg(end) - Time_years_reg(1);
            total_mass_change_reg = total_mass_rate_reg * time_span_reg;

            fprintf('\n==============================\n');
            fprintf('Regression-based mass check: %s\n', data_name);
            fprintf('Time span: %.2f years\n', time_span_reg);
            fprintf('Mass rate: %.1f Gt/yr\n', total_mass_rate_reg);
            fprintf('Total mass change: %.1f Gt\n', total_mass_change_reg);
            fprintf('==============================\n');
        end
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
        fname_elas = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/altimetry/khan_et_al_2025/Greenland_dhdt_elas_1kmgrid_DB.nc';
        fname_GIA = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/altimetry/khan_et_al_2025/Greenland_GIA_1kmgrid_DB.nc';

        % Read the variables
        X = ncread(filename, 'x');       % East coordinate (km)
        Y = ncread(filename, 'y');       % North coordinate (km)
        Time = ncread(filename, 'time'); % Time in days since 2003-01-01
        dhdt_vol = ncread(filename, 'dhdt_vol'); % Mean monthly elevation change rate (m/month) in ice equivalent height, fill value is -9999
        dhdt_firn = ncread(fname_firn, 'dhdt_firn'); % Mean monthly firn compaction rate (m/month), fill value is -9999
        dhdt_elas = ncread(fname_elas, 'dhdt_elas'); % Elastic uplift correction in mm/month (time, y, x), fill value is -9999
        dhdt_GIA = ncread(fname_GIA, 'dhdt_GIA'); % GIA correction in mm/yr (y, x), fill value is NaN

        % Convert fill values to NaN
        dhdt_vol(dhdt_vol == -9999) = NaN;
        dhdt_firn(dhdt_firn == -9999) = NaN;
        dhdt_elas(dhdt_elas == -9999) = NaN;

        % Convert units in mm/time to m/month
        dhdt_elas = dhdt_elas * 1e-3;
        dhdt_GIA = (1/12) * dhdt_GIA * 1e-3;

        % Calculate the total monthly change in ice elevation (m/month)
        dhdt_monthly = dhdt_vol - dhdt_firn + dhdt_elas + dhdt_GIA;

        % Convert fill values to NaN (redundant but to be safe)
        dhdt_monthly(dhdt_monthly == -9999) = NaN;

        % Calculate monthly ice mass change
        dhdt_mass_monthly = 0.917 * dhdt_vol ...
                          - dhdt_firn ...
                          + dhdt_elas ...
                          + dhdt_GIA;

        % Convert time from days since 2003-01-01 to decimal years
        reference_date = datetime(2003, 1, 1);
        Time_datetime = reference_date + days(Time);
        Time = year(Time_datetime) + (day(Time_datetime, 'dayofyear') - 1) / 365.25;

        % Filter data to end of 2022 (exclude incomplete 2023 data)
        time_mask = Time <= 2023.0; % Include only data up to end of 2022
        Time = Time(time_mask);
        dhdt_monthly = dhdt_monthly(:, :, time_mask);
        dhdt_mass_monthly = dhdt_mass_monthly(:,:, time_mask);

        %% Note: If want to match Khan et al.'s calculation, do not mask out the dhdt_monthly field, and run
        % [dice_mass, dice_mass_cum] = get_ice_mass_change(1000, X, Y, dhdt_mass_monthly(:,:,2:end), Time(2:end), false);
        % This way, it matches their calculation:
        %     dhdt_mass_monthly(isnan(dhdt_mass_monthly))=0;
        %     acc_mass(1)=0;
        %     for i=2:247
        %        acc_mass(i)=acc_mass(i-1)+0.001*sum(sum(dhdt_mass_monthly(:,:,i)));
        %     end
        % The answer from this calculation is -4538.8 Gt.

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
    if strcmp(data_name, 'measureItsLive-GEMB') || strcmp(data_name, 'measureItsLive-GSFC') || strcmp(data_name, 'DTU2016')
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
        dhdt_annual_raw = squeeze(sum(dhdt_monthly_reshaped, 3, 'omitnan')); % Before quality check
        
        % Sum over the 12 months (3rd dimension) to get annual changes ignoring NaNs
        % Require enough valid months per year to avoid totals dominated by incomplete temporal coverage.
        % That is, for each pixel and each individual year, fewer than 10 out of 12 monthly values are valid,
        % the annual value for that pixel-year is discarded (set to NaN).
        valid_count = sum(~isnan(dhdt_monthly_reshaped), 3);

        dhdt_annual_qc = sum(dhdt_monthly_reshaped, 3, 'omitnan');
        dhdt_annual_qc(valid_count < 10) = NaN; % require at least 10 valid months

        dhdt_annual_qc = squeeze(dhdt_annual_qc);
        dhdt_annual = dhdt_annual_qc;
    elseif strcmp(data_name, 'Buffalo2025-GEMB') || strcmp(data_name, 'Buffalo2025-GSFC') || strcmp(data_name, 'Buffalo2025-IMAU')
        num_years = size(dhdt_annual, 3);
        % Other than that,do nothing since the data is already in annual intervals
        dhdt_monthly = NaN;
    end
    
    %% DEBUG: inspect dhdt_annual on ORIGINAL altimetry grid
    fprintf('\n=====================================\n');
    fprintf('DEBUG ORIGINAL ALT GRID: %s\n', data_name);
    fprintf('Size X: %s\n', mat2str(size(X)));
    fprintf('Size Y: %s\n', mat2str(size(Y)));
    fprintf('Size dhdt_annual: %s\n', mat2str(size(dhdt_annual)));
    fprintf('Years: %.2f to %.2f, n=%d\n', years(1), years(end), numel(years));

    tmp = dhdt_annual;
    fprintf('NaNs: %d / %d = %.2f%%\n', ...
        sum(isnan(tmp(:))), numel(tmp), 100*sum(isnan(tmp(:)))/numel(tmp));

    fprintf('Range original dhdt_annual: %.3f to %.3f m/yr\n', ...
        min(tmp(:), [], 'omitnan'), max(tmp(:), [], 'omitnan'));

    fprintf('Mean/std original dhdt_annual: %.3f / %.3f m/yr\n', ...
        mean(tmp(:), 'omitnan'), std(tmp(:), 'omitnan'));

    % Per-year range check
    for kk = 1:size(dhdt_annual,3)
        d = dhdt_annual(:,:,kk);
        fprintf('Year %.0f: min %.3f, max %.3f, mean %.3f, std %.3f m/yr\n', ...
            years(kk), ...
            min(d(:), [], 'omitnan'), ...
            max(d(:), [], 'omitnan'), ...
            mean(d(:), 'omitnan'), ...
            std(d(:), 'omitnan'));
    end

    % Find extreme pixels
    %abs_thresh = 20; % m/yr; adjust if needed
    %bad = abs(dhdt_annual) > abs_thresh;
    %fprintf('Pixels with |dhdt_annual| > %.1f m/yr: %d\n', ...
    %    abs_thresh, sum(bad(:), 'omitnan'));

    %[i_bad, j_bad, t_bad] = ind2sub(size(dhdt_annual), find(bad, 10, 'first'));

    %for nn = 1:numel(i_bad)
    %    fprintf('Bad pixel #%d: i=%d, j=%d, year=%.0f, dhdt=%.3f m/yr\n', ...
    %        nn, i_bad(nn), j_bad(nn), years(t_bad(nn)), ...
    %        dhdt_annual(i_bad(nn), j_bad(nn), t_bad(nn)));
    %end

    % Plot one suspicious year on original X-Y grid
    %[~, worst_t] = max(squeeze(max(max(abs(dhdt_annual), [], 1, 'omitnan'), [], 2, 'omitnan')));

    %figure('Color','w');
    %imagesc(X, Y, dhdt_annual(:,:,worst_t)');
    %set(gca,'YDir','normal');
    %axis equal tight;
    %colorbar;
    %title(sprintf('%s original grid dhdt annual, year %.0f', data_name, years(worst_t)));
    %xlabel('X');
    %ylabel('Y');
    %caxis([-5 5]); % use wider range if needed
    fprintf('=====================================\n\n');

    % Calculate the total change across all years, ignoring NaNs
    dhdt_total = sum(dhdt_annual, 3, 'omitnan');

    % Convert years to double precision for polyfit compatibility in get_ice_mass_change function
    years = double(years);

    % Create years_thickness for thickness data (has nt+1 elements)
    switch data_name
        case {'measureItsLive-GEMB', 'measureItsLive-GSFC', 'DTU2016', 'DTU2025'}
            % These datasets end mid- or end-year; extend by +1 year
            years_thickness = [years; years(end) + 1];

        case {'Buffalo2025-GEMB', 'Buffalo2025-GSFC', 'Buffalo2025-IMAU'}
            years_thickness = [years(1) - 1; years];
end

% Calculate ice mass change in Gt
if strcmp(data_name, 'DTU2025')
    disp("Calculating ice mass change based on monthly mass-equivalent field")
    [dice_mass_monthly_mass, dice_mass_cum_mass] = get_ice_mass_change(rhoo, X, Y, dhdt_mass_monthly, Time, false);
    disp(" ")

    unique_years_mass = unique(floor(Time));
    num_years_mass = length(unique_years_mass);
    dhdt_mass_annual = zeros(size(dhdt_mass_monthly, 1), size(dhdt_mass_monthly, 2), num_years_mass);

    for i = 1:num_years_mass
        year_mask = floor(Time) == unique_years_mass(i);
        if any(year_mask)
            dhdt_mass_annual(:,:,i) = sum(dhdt_mass_monthly(:,:,year_mask), 3, 'omitnan');
        end
    end

    dhdt_mass_total = sum(dhdt_mass_annual, 3, 'omitnan');

    disp("Calculating ice mass change based on annualized mass-equivalent field")
    [dice_mass_annual_mass, dice_mass_cum_annual_mass] = get_ice_mass_change(rhoo, X, Y, dhdt_mass_annual, unique_years_mass, false);
    disp(" ")

    disp("Calculating ice mass change based on total mass-equivalent field")
    [dice_mass_total_mass, dice_mass_cum_total_mass] = get_ice_mass_change(rhoo, X, Y, dhdt_mass_total, 1, false);
    disp(" ")

elseif strcmp(data_name, 'Buffalo2025-GEMB') || strcmp(data_name, 'Buffalo2025-GSFC') || strcmp(data_name, 'Buffalo2025-IMAU')
    rho = rhoi;
    disp("Calculating ice mass change based on annual elevation-change data")
    [dice_mass_annual, dice_mass_cum] = get_ice_mass_change(rho, X, Y, dhdt_annual, years, false);
    disp(" ")

    disp("Calculating ice mass change based on the total elevation change across data coverage")
    [dice_mass_total, dice_mass_cum_total] = get_ice_mass_change(rho, X, Y, dhdt_total, 1, false);
    disp(" ")
elseif strcmp(data_name, 'measureItsLive-GEMB') || strcmp(data_name, 'measureItsLive-GSFC')
    rho = rhoi;
    disp("Calculating ice mass change based on monthly elevation-change field")
    [dice_mass_monthly_mass, dice_mass_cum_mass] = get_ice_mass_change(rho, X, Y, dhdt_monthly, Time, false);
    disp(" ")

    disp("Calculating ice mass change based on annual elevation-change data - BEFORE removing year with less than 10 valid months")
    [dice_mass_annual, dice_mass_cum] = get_ice_mass_change(rho, X, Y, dhdt_annual_raw, years, false);
    disp(" ")

    disp("Calculating ice mass change based on annual elevation-change data - AFTER removing year with less than 10 valid months")
    [dice_mass_annual, dice_mass_cum] = get_ice_mass_change(rho, X, Y, dhdt_annual, years, false);
    disp(" ")

    disp("Calculating ice mass change based on the total elevation change across data coverage - After QC")
    [dice_mass_total, dice_mass_cum_total] = get_ice_mass_change(rho, X, Y, dhdt_total, 1, false);
    disp(" ")
end


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
%[lat_sphere, long_sphere, dhdt_annual] = ellipsoid_to_sphere(lat_ellipsoid, long_ellipsoid, r_earth, dhdt_annual);

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
fprintf('Std dhdt_annual: %.3f m/yr\n', std(dhdt_annual(:), 'omitnan'));

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

if strcmp(data_name, 'measureItsLive-GEMB') || strcmp(data_name, 'measureItsLive-GSFC') || strcmp(data_name, 'DTU2016')
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
