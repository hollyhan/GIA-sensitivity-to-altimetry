function [lat_matrix, lon_matrix, uu_detrended_cell, uu_err, time_gnss, R2_gnss, rate_gnss] = preprocess_gnss_data(stn_id, fpath_gnss, fname_coord_gnss, n_degree)
% PREPROCESS_GNSS_DATA Processes GNSS data for multiple stations and returns detrended vertical signal with coordinates.
%   Inputs:
%       stn_id            - cell array of station IDs (e.g., {'KELY','QAQ1'})
%       fpath_gnss        - path to folder containing GNSS .txt files
%       fname_coord_gnss  - filename of coordinate file
%       n_degree          - degree of polynomial for detrending
%   Outputs:
%       lat_matrix        - matrix of latitudes for all stations
%       lon_matrix        - matrix of longitudes for all stations
%       uu_detrended_cell - cell array of detrended vertical signals for all stations
%       uu_err            - cell array or gnss error for all stations
%       time_gnss         - cell array of time stamps for all stations
%       R2_gnss           - vector array of R-squared values for all stations
%       rate_gnss         - vector array of rates of change for all stations

    % Load settings from an external file
    run('settings_observation_data.m')

    % Read GNSS filenames and coordinates file
    fnames_gnss = dir(fullfile(fpath_gnss, '*.txt'));
    coord_data_gnss = readlines(fullfile(fpath_gnss, fname_coord_gnss));
    coord_data_gnss = coord_data_gnss(~startsWith(coord_data_gnss, '%')); % Remove comments

    % Initialize matrices for latitudes, longitudes, and a cell array for detrended vertical data
    num_stations = length(stn_id);
    lat_matrix = NaN(num_stations, 1);  % Matrix for latitudes of all stations
    lon_matrix = NaN(num_stations, 1);  % Matrix for longitudes of all stations
    uu_detrended_cell = cell(num_stations, 1);  % Cell array for detrended vertical data
    uu_err = cell(num_stations, 1); % Cell array for error at each site
    time_gnss = cell(num_stations, 1);
    rate_gnss = NaN(num_stations, 1);  % Initialize rate output

    % Loop over each station
    my_color = jet(num_stations);  % Color map for each station
    for ii = 1:num_stations

        % Prepare figure for plotting
        figure('Position', [600 1500 550 250])
        set(0, 'DefaultAxesFontSize', 17, ...
            'DefaultAxesLineWidth', 1, ...
            'DefaultTextFontSize', 17, ...
            'DefaultLineMarkerSize', 8);

        % Plot settings
        axes1 = axes('Layer', 'top', 'Position', [0.1 0.1 0.8 0.8]);
        box(axes1, 'on'); hold(axes1, 'all'); grid on;
        ylabel(axes1, 'Up-down motion [mm]');

        % Find the corresponding GNSS data file for the station
        for kk = 1:length(fnames_gnss)
            name = fnames_gnss(kk).name;
            if strcmp(stn_id{ii}, name(1:4))
                % Get coordinates from the coordinate data file
                lat = NaN; lon = NaN;
                for i = 1:length(coord_data_gnss)
                    parts = strsplit(strtrim(coord_data_gnss(i)));
                    if length(parts) >= 3 && strcmp(stn_id{ii}, parts{1})
                        lat = str2double(parts{2});
                        lon = str2double(parts{3});
                        break;
                    end
                end

                % Store latitude and longitude in the matrix
                lat_matrix(ii) = lat;
                lon_matrix(ii) = lon;

                % Load GNSS data for the station
                filepath = fullfile(fpath_gnss, name);
                disp(['Station to analyze: ', name(1:4)]);

                % Read the data, assuming first 4 lines are headers
                data = readmatrix(filepath, 'NumHeaderLines', 4);
                tt = data(:, 1);   % Time [year]
                uu = data(:, 6);   % Vertical displacement [mm]
                uu_err_tmp = data(:, 7); % Error in vertical displacement [mm]

                % Get weighted annual average vertical displacement and uncertainty
                if annual_output
                    unique_yr = unique(floor(tt));

                     % Initialize arrays
                    annual_avg = NaN(length(unique_yr), 1);
                    annual_err_mean = NaN(length(unique_yr), 1);
                    annual_err = NaN(length(unique_yr), 1);
                    annual_std = NaN(length(unique_yr), 1);

                    for yy = 1:length(unique_yr)
                        yr = unique_yr(yy);
                        idx_unique_yr = find(floor(tt) == yr);

                        if ~isempty(idx_unique_yr)
                            val = uu(idx_unique_yr);
                            sigma = uu_err_tmp(idx_unique_yr);
                            w = 1 ./ sigma.^2;

                            annual_avg(yy) = sum(w .* val) / sum(w); % weighted mean
                            annual_err_mean(yy) = sqrt(1 / sum(w)); % Formal uncertainty of weighted mean
                            annual_std(yy) = std(val, 'omitnan'); % standard deviation of the data. More appropriate than the range of the data for weighted mean
                            annual_err(yy) = max(annual_err_mean(yy), annual_std(yy)); % pick the larger of the formal error and the empirical standard deviation
                        end
                    end

                    % Detrend the vertical displacement data
                    fprintf('Detrending the GNSS annual averages with polynomial degree %d\n', n_degree);
                    if n_degree == 0
                        disp('Removing mean from the signal');
                        uu_detrended = detrend(annual_avg, 0);
                    else
                        uu_detrended = detrend(annual_avg, n_degree);
                    end
                    time_gnss{ii} = unique_yr;
                    uu_err{ii} = annual_err;

                    % Fit linear trend to get rate of change
                    mdl_linear = fitlm(unique_yr, uu_detrended);
                    rate_gnss(ii) = mdl_linear.Coefficients.Estimate(2);      % Slope (mm/yr)
                    intercept = mdl_linear.Coefficients.Estimate(1); % Intercept
                    y_fit = mdl_linear.Fitted;                       % Fitted values
                    R2 = mdl_linear.Rsquared.Ordinary;
                    R2_gnss(ii) = R2;
                else
                    % Detrend the vertical displacement data
                    fprintf('Detrending the GNSS data with polynomial degree %d\n', n_degree);
                    if n_degree == 0
                        disp('Removing mean from the signal');
                        uu_detrended = detrend(uu, 0);
                    else
                        uu_detrended = detrend(uu, n_degree);
                    end
                    time_gnss{ii} = tt;
                    uu_err{ii} = uu_err_tmp;

                    % Fit linear trend to get rate of change
                    mdl_linear = fitlm(tt, uu_detrended);
                    rate_gnss(ii) = mdl_linear.Coefficients.Estimate(2);      % Slope (mm/yr)
                    intercept = mdl_linear.Coefficients.Estimate(1); % Intercept
                    y_fit = mdl_linear.Fitted;                       % Fitted values
                    R2 = mdl_linear.Rsquared.Ordinary;
                    R2_gnss(ii) = R2;
                end

                % Store detrended vertical data in the cell array
                uu_detrended_cell{ii} = uu_detrended;

                % Test a degree of linearity in the signal
                if annual_output
                    % Shade the uncertainty region
                    fill([unique_yr; flipud(unique_yr)], ...
                        [uu_detrended - uu_err{ii}; flipud(uu_detrended + uu_err{ii})], ...
                        my_color(ii,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'Parent', axes1);
                    hold on
                    % Then plot the data points
                    plot(unique_yr, uu_detrended, 'b.', 'Parent', axes1, 'MarkerSize',10);
                    fit_label = sprintf('Linear fit (%.2f mm/yr)', rate_gnss(ii));
                    plot(unique_yr, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', fit_label);
                    legend("Detrended Average Data", fit_label, 'Location', 'southeast');
                else
                    % Shade the uncertainty region
                    %fill([tt; flipud(tt)], ...
                    %    [uu_detrended - uu_err_tmp; flipud(uu_detrended + uu_err_tmp)], ...
                    %    my_color(ii,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'Parent', axes1);
                    hold on
                    % Then plot the data points
                    plot(tt, uu_detrended, '.r');
                end

                title(['Station: ', stn_id{ii}(1:4), ...
                       '   |   Linearity (R²): ', num2str(R2, '%.2f')]);
                set(gcf, 'color', 'w');

                if savefig_gnss
                    if ~exist('./Figures', 'dir'); mkdir('./Figures'); end
                    filename = fullfile('./Figures', ['GNET_' stn_id{ii} '_vertical.png']);
                    exportgraphics(gcf, filename, 'Resolution', 300);
                end
            end
        end
    end
end