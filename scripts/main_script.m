%% Define which steps to run
steps=[1];

addpath('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/scripts/mesh');
% Load settings
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
    [h_annual_1, dhdt_annual_1, dhdt_monthly_1, years_altimetry_1, lat_sphere_1, long_sphere_1, X_1, Y_1] = preprocess_ice_altimetry('measureItsLive-GEMB', false);
    [h_annual_2, dhdt_annual_2, dhdt_monthly_2, years_altimetry_2, lat_sphere_2, long_sphere_2, X_2, Y_2] = preprocess_ice_altimetry('measureItsLive-GSFC', false);
    [h_annual_3, dhdt_annual_3, dhdt_monthly_3, years_altimetry_3, lat_sphere_3, long_sphere_3, X_3, Y_3] = preprocess_ice_altimetry('DTU2016', false);
    [h_annual_4, dhdt_annual_4, dhdt_monthly_4, years_altimetry_4, lat_sphere_4, long_sphere_4, X_4, Y_4] = preprocess_ice_altimetry('DTU2025', false);% DTU data reports-4186.2778 Gt between 2003-2022-12-31 and 4701 Gt if not correcting for firn, Elastic uplift and GIA
    [h_annual_5, dhdt_annual_5, dhdt_monthly_5, years_altimetry_5, lat_sphere_5, long_sphere_5, X_5, Y_5] = preprocess_ice_altimetry('Buffalo2025-GEMB', false);
    [h_annual_6, dhdt_annual_6, dhdt_monthly_6, years_altimetry_6, lat_sphere_6, long_sphere_6, X_6, Y_6] = preprocess_ice_altimetry('Buffalo2025-GSFC', false);
    [h_annual_7, dhdt_annual_7, dhdt_monthly_7, years_altimetry_7, lat_sphere_7, long_sphere_7, X_7, Y_7] = preprocess_ice_altimetry('Buffalo2025-IMAU', false);
    disp('====================================');

    years_altimetry = [years_altimetry_1; years_altimetry_2; years_altimetry_3; years_altimetry_4; years_altimetry_5; years_altimetry_6; years_altimetry_7];
    years_altimetry = unique(years_altimetry);

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
    % Process firn model datasets
    disp('=== Processing firn model data ===');
    [dfac_annual_1, dfac_monthly_1, years_firn_1, lat_sphere_firn_1, long_sphere_firn_1] = preprocess_firn_model('GEMB'); % preprocessed for measureItsLive
    [dfac_annual_2, dfac_monthly_2, years_firn_2, lat_sphere_firn_2, long_sphere_firn_2] = preprocess_firn_model('GSFC'); % preprocessed for measureItsLive
    [dfac_annual_3, dfac_monthly_3, years_firn_3, lat_sphere_firn_3, long_sphere_firn_3] = preprocess_firn_model('RACMO2.3p2'); % preprocessed for DTU2025
    disp('====================================');
end

if any(steps==4)
    % Process glacier mask datasets (returns native high-resolution mask on spherical geographic coordinates)
    disp('=== Processing glacier mask data and apply the mask to the altimetry datasets ===');
    disp('---- Altimetry Dataset: MeaSURES ITS_LIVE ----')

    % Find overlapping years between the mask and altimetry data
    yrs_mask = 1985:2022; % manually defined reliable data years based on the metadata
    yrs_overlap = intersect(years_altimetry_1, yrs_mask);
    disp(['Dataset 1 (measureItsLive): mask years ', num2str(yrs_overlap(1)), '-', num2str(yrs_overlap(end))]);

    % Process mask year-by-year to avoid memory issues
    h_masked = zeros(size(h_annual_1,1), size(h_annual_1,2), length(yrs_overlap));
    for n = 1:length(yrs_overlap)
        % Extract the native mask data
        [mask, mask_year(n), lat_mask, lon_mask, x_mask, y_mask] = preprocess_glacier_mask(yrs_overlap(n));

        % Resample onto the altimetry grid
        disp('=== Resampling glacier mask to each altimetry dataset ==='); % This part needs to be debugged because ice_masks are processed in parts
        mask_resampled = resample_mask_to_target_grid(mask, lat_mask, lon_mask, lat_sphere_1, long_sphere_1);
        mask_resampled = double(mask_resampled); % convert from logical to numeric array

        % Mask thickness data
        h_masked(:,:,n) = h_annual_1(:,:,n).*mask_resampled;
    end

    % Derive dhdt from h_masked
    dhdt_masked = zeros(size(h_masked,1), size(h_masked,2), length(yrs_overlap)-1);
    for m = 1:length(yrs_overlap)-1
        dhdt_masked(:,:,m) = h_masked(:,:,m+1) - h_masked(:,:,m); % This is the same as dhdt_annual_1 as benchmarked in step 2

        plot_mask = false;
        if plot_mask
            np = 64; % Number of colors
            blue_to_white = [linspace(0,1,np/2)', linspace(0,1,np/2)', ones(np/2,1)];
            white_to_red = [ones(np/2,1), linspace(1,0,np/2)', linspace(1,0,np/2)'];
            custom_cmap = [blue_to_white; white_to_red];

            % For pcolor, we need to extract 1D coordinate vectors
            % Take the first row and first column as representative vectors
            lat_vec = lat_sphere_1(:,1);  % First column as latitude vector
            lon_vec = long_sphere_1(1,:); % First row as longitude vector

            figure;
            data = dhdt_masked(:,:,n);
            pcolor(lon_vec,lat_vec, data);
            %pcolor(lon_vec, lat_vec, mask_resampled);
            shading flat;
            colorbar;
            title(sprintf('masked dhdt between year %d and %d', yrs_overlap(m), yrs_overlap(m+1)));
            xlabel('Longitude', 'FontSize', 14);
            ylabel('Latitude', 'FontSize', 14);
            set(gca, 'FontSize', 14);
            colormap(flip(custom_cmap))%
            caxis([-1 1])
        end
    end
    disp('====================================');
end


if any(steps==5)
    addpath('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/scripts/mesh');
    % Load mesh model and initialize models for each altimetry dataset
    md_regional = loadmodel(fpath_mesh_model_regional);  % can be anything (e.g., '[]') if bool_mesh_greenland_external is 'false'

    % Interpolate with mass conservation enabled and corresponding ice masks
    [md1_regional, mass_report_1] = interpolate_altimetry_to_mesh(h_annual_1, lat_sphere_1, long_sphere_1, years_altimetry_1, md_regional, X_1, Y_1, dhdt_annual_1);
    [md1, ~, ~, ~] = createGlobalMesh(md1_regional); % create global mesh from regional mesh
    md1 = initialize_model(md1);

    [md2_regional, mass_report_2] = interpolate_altimetry_to_mesh(h_annual_2, lat_sphere_2, long_sphere_2, years_altimetry_2, md_regional, X_2, Y_2, dhdt_annual_2);
    [md2, ~, ~, ~] = createGlobalMesh(md2_regional); % create global mesh from regional mesh
    md2 = initialize_model(md2);

    [md3_regional, mass_report_3] = interpolate_altimetry_to_mesh(h_annual_3, lat_sphere_3, long_sphere_3, years_altimetry_3, md_regional, X_3, Y_3, dhdt_annual_3);
    [md3, ~, ~, ~] = createGlobalMesh(md3_regional); % create global mesh from regional mesh
    md3 = initialize_model(md3);

    [md4_regional, mass_report_4] = interpolate_altimetry_to_mesh(h_annual_4, lat_sphere_4, long_sphere_4, years_altimetry_4, md_regional, X_4, Y_4, dhdt_annual_4);
    [md4, ~, ~, ~] = createGlobalMesh(md4_regional); % create global mesh from regional mesh
    md4 = initialize_model(md4);

    [md5_regional, mass_report_5] = interpolate_altimetry_to_mesh(h_annual_5, lat_sphere_5, long_sphere_5, years_altimetry_5, md_regional, X_5, Y_5, dhdt_annual_5);
    [md5, ~, ~, ~] = createGlobalMesh(md5_regional); % create global mesh from regional mesh
    md5 = initialize_model(md5);

    [md6_regional, mass_report_6] = interpolate_altimetry_to_mesh(h_annual_6, lat_sphere_6, long_sphere_6, years_altimetry_6, md_regional, X_6, Y_6, dhdt_annual_6);
    [md6, ~, ~, ~] = createGlobalMesh(md6_regional); % create global mesh from regional mesh
    md6 = initialize_model(md6);

    [md7_regional, mass_report_7] = interpolate_altimetry_to_mesh(h_annual_7, lat_sphere_7, long_sphere_7, years_altimetry_7, md_regional, X_7, Y_7, dhdt_annual_7);
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

    figure()
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
    %md1_solved = run_gia(md1, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md1_solved');

    % with ice profile #2
    %[md2, vlm2_VE, vlm2_elastic, hlm2, accm2] = run_gia_greensFunction(md2, ht, lt, kt, false, lat_gnss, lon_gnss);
    %md2_solved = run_gia(md2, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md2_solved');

    % with ice profile #3
    %[md3, vlm3_VE, vlm3_elastic, hlm3, accm3] = run_gia_greensFunction(md3, ht, lt, kt, false, lat_gnss, lon_gnss);
    %md3_solved = run_gia(md3, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md3_solved');

    % with ice profile #4
    %[md4, vlm4_VE, vlm4_elastic, hlm4, accm4] = run_gia_greensFunction(md4, ht, lt, kt, false, lat_gnss, lon_gnss);
    md4_solved = run_gia(md4, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md4_solved');

    % with ice profile #5
    %[md5, vlm5_VE, vlm5_elastic, hlm5, accm5] = run_gia_greensFunction(md5, ht, lt, kt, false, lat_gnss, lon_gnss);
    %md5_solved = run_gia(md5, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md5_solved');

    % with ice profile #6
    %[md6, vlm6_VE, vlm6_elastic, hlm6, accm6] = run_gia_greensFunction(md6, ht, lt, kt, false, lat_gnss, lon_gnss);
    md6_solved = run_gia(md6, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md6_solved');

    % with ice profile #7
    %[md7, vlm7_VE, vlm7_elastic, hlm7, accm7] = run_gia_greensFunction(md7, ht, lt, kt, false, lat_gnss, lon_gnss);
    md7_solved = run_gia(md7, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, 'md7_solved');

    disp('====================================');
end



addpath('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results/model_objects_saved/')
md1_solved = loadmodel('md1_solved.mat');
md2_solved = loadmodel('md2_solved.mat');
md3_solved = loadmodel('md3_solved.mat');
md5_solved = loadmodel('md5_solved.mat');
md6_solved = loadmodel('md6_solved.mat');
md7_solved = loadmodel('md7_solved.mat');

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
    saveas(gcf, sprintf('residual_comparison_for_different_ice_loading_models-elastic.png'));
end