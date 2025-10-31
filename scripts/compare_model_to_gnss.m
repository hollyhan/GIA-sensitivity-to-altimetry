function [WRMS, mean_err_gnss, interped_md_data, misfit, rate, rate_gnss_fit, y_fit_model_all, y_fit_gnss_all] = compare_model_to_gnss(lat_gnss, lon_gnss, data_gnss, err_gnss, common_time, stn_id, md, vlm_GF_gnss)
%   Inputs:
%       lat_gnss            - vector array of latitudes for all stations
%       lon_gnss            - vector array of longitudes for all stations
%       data_gnss           - cell array of detrended vertical GNSS signals for all stations
%       err_gnss            - cell array or gnss error for all stations
%       common_time         - cell array of time stamps for model-data comparison
%       stn_id              - cell array of station IDs for all stations
%       md                  - model object
%       vlm_GF_gnss         - vertical landmotion at gnss sites, if applicable (use_whole_domain = false)
%
%   Outputs:
%       WRMS                - weighted root mean square error between the model and gnss data
%       mean_err_gnss       - vector array of mean error per GNSS station
%       interped_md_data    - interpolated (spatiotemporally) data_field from the input model data_field
%       misfit              - vector array of misfit between the model and gnss data
%       rate                - vector array of rate of change of the model
%       rate_gnss_fit       - vector array of rate of change of the gnss data
%       y_fit_model_all     - cell array of fitted model data
%       y_fit_gnss_all      - cell array of fitted gnss data

    % Check if we're using whole domain or GNSS sites only
    run('settings_gia_parameterization.m')
    plot_fig = false;

    % Initialize arrays for rates and fit data
    rate = zeros(length(data_gnss), 1);
    rate_gnss_fit = zeros(length(data_gnss), 1);
    y_fit_model_all = cell(length(data_gnss), 1);
    y_fit_gnss_all = cell(length(data_gnss), 1);

    % Interpolate data onto GNSS sites and common timesteps, and detrend
    if ~isempty(md)
        % Get model outputs
        [~, ~, ~, ~, delta_B, ~, time_md_pres] = get_GIAoutputs(md);

        if use_whole_domain
            % Need to interpolate from whole domain
            disp('Interpolating present-day model data field onto GPS sites')
            interped_matrix = interp_model_to_gnss(md, delta_B, lat_gnss, lon_gnss, 'nearest');

            for n = 1:length(data_gnss)
                % Interpolate onto common times
                interped_md_data{n} = interp1(time_md_pres, interped_matrix(n,:), common_time{n}, 'linear', 'extrap');
                %interped_md_data{n} = interped_matrix(n,:); % use this if you ever want to use the original model data
                % Remove mean (no need to detrend linearly because no GIA signal is included)
                disp('Removing mean from the interpolated data')
                interped_md_data{n} = detrend(interped_md_data{n}, 0);

                % Convert to mm, consistent with the GNSS data
                interped_md_data{n} = interped_md_data{n} * 1e3;
                % Fit linear trend to get rate of change using fitlm
                mdl_linear = fitlm(common_time{n}, interped_md_data{n});
                rate(n) = mdl_linear.Coefficients.Estimate(2);  % Rate of change in mm/yr (slope)
                intercept{n} = mdl_linear.Coefficients.Estimate(1);  % Intercept in mm

                % Calculate GNSS data linear fit (always needed for misfit calculation)
                mdl_gnss = fitlm(common_time{n}, data_gnss{n});
                rate_gnss_fit(n) = mdl_gnss.Coefficients.Estimate(2);  % Rate of change in mm/yr

                % Store fit data for output
                y_fit_model_all{n} = mdl_linear.Fitted;
                y_fit_gnss_all{n} = mdl_gnss.Fitted;

                if plot_fig
                    % Create figure showing both model and GNSS data with slope lines
                    figure;

                    % Plot model data and its linear fit
                    plot(common_time{n}, interped_md_data{n}, 'b.', 'DisplayName', 'Model Data', 'MarkerSize', 8);
                    hold on;
                    y_fit_model = mdl_linear.Fitted;
                    plot(common_time{n}, y_fit_model, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Model fit (%.2f mm/yr)', rate{n}));

                    % Calculate and plot GNSS data linear fit
                    y_fit_gnss = mdl_gnss.Fitted;
                    R2_gnss = mdl_gnss.Rsquared.Ordinary;

                    % Plot GNSS data and its linear fit
                    plot(common_time{n}, data_gnss{n}, 'r.', 'DisplayName', 'GNSS Data', 'MarkerSize', 8);
                    plot(common_time{n}, y_fit_gnss, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('GNSS fit (%.2f mm/yr)', rate_gnss_fit{n}));

                    xlabel('Time (year)');
                    ylabel('VLM (mm)');
                    title(sprintf('VLM Comparison at Station %s (Model R^2 = %.3f, GNSS R^2 = %.3f)', stn_id{n}, mdl_linear.Rsquared.Ordinary, R2_gnss));
                    legend('show', 'Location', 'best');
                    grid on;
                    set(gca, 'FontSize', 12);
                end
            end
        else
            % Data is already at GNSS sites
            disp('Using model results at GNSS sites')
            for n = 1:length(data_gnss)
                % Interpolate onto common times if needed - maybe we don't need this
                if ~isequal(time_md_pres, common_time{n})
                    interped_md_data{n} = interp1(time_md_pres, vlm_GF_gnss(n,:), common_time{n}, 'linear', 'extrap');
                else
                    interped_md_data{n} = vlm_GF_gnss(n,:);
                end

                % Remove mean
                disp('Removing mean from the interpolated data')
                interped_md_data{n} = detrend(interped_md_data{n}, 0);

                % Convert to mm
                interped_md_data{n} = interped_md_data{n} * 1e3;

                % Fit linear trend to get rate of change using fitlm
                mdl_linear = fitlm(common_time{n}, interped_md_data{n});
                rate(n) = mdl_linear.Coefficients.Estimate(2);  % Rate of change in mm/yr (slope)
                intercept{n} = mdl_linear.Coefficients.Estimate(1);  % Intercept in mm

                % Calculate GNSS data linear fit (always needed for misfit calculation)
                mdl_gnss = fitlm(common_time{n}, data_gnss{n});
                rate_gnss_fit(n) = mdl_gnss.Coefficients.Estimate(2);  % Rate of change in mm/yr

                % Store fit data for output
                y_fit_model_all{n} = mdl_linear.Fitted;
                y_fit_gnss_all{n} = mdl_gnss.Fitted;

                % Get R-squared from the model
                R2 = mdl_linear.Rsquared.Ordinary;

                if plot_fig
                    % Create figure showing both model and GNSS data with slope lines
                    figure;

                    % Plot model data and its linear fit
                    plot(common_time{n}, interped_md_data{n}, 'b.', 'DisplayName', 'Model Data', 'MarkerSize', 8);
                    hold on;
                    y_fit_model = mdl_linear.Fitted;
                    plot(common_time{n}, y_fit_model, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Model fit (%.2f mm/yr)', rate{n}));

                    % Calculate and plot GNSS data linear fit
                    y_fit_gnss = mdl_gnss.Fitted;
                    R2_gnss = mdl_gnss.Rsquared.Ordinary;

                    % Plot GNSS data and its linear fit
                    plot(common_time{n}, data_gnss{n}, 'r.', 'DisplayName', 'GNSS Data', 'MarkerSize', 8);
                    plot(common_time{n}, y_fit_gnss, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('GNSS fit (%.2f mm/yr)', rate_gnss_fit{n}));

                    xlabel('Time (year)');
                    ylabel('VLM (mm)');
                    title(sprintf('VLM Comparison at Station %s (Model R^2 = %.3f, GNSS R^2 = %.3f)', stn_id{n}, mdl_linear.Rsquared.Ordinary, R2_gnss));
                    legend('show', 'Location', 'best');
                    grid on;
                    set(gca, 'FontSize', 12);
                end
            end
        end
    end

    % Calculate WRMS
    [WRMS, mean_err_gnss] = calculate_wrms_and_meanerr(data_gnss, err_gnss, interped_md_data);

    % calculate misfit for each station with different time series
    misfit = zeros(length(data_gnss), 1);
    for n = 1:length(data_gnss)
        % Compute normalized residuals for station n
        misfit(n)= (rate(n) - rate_gnss_fit(n)) ./ mean(err_gnss{n});
    end
end
