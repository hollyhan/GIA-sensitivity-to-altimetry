% load md_solved for ISSM elastic simulation
md_ISSM = loadmodel('/Users/kyhan/Desktop/Projects/Greenland-LIA/results/model_objects_saved/md_gia_solved_present_1993-2023CE_Elastic_iter0.mat');
md_ISSM = loadmodel('/Users/kyhan/Desktop/Projects/Greenland-LIA/results/model_objects_saved/md_gia_solved_present_1993-2023CE_Maxwell_iter0.mat'); %for Viscoelastic case

% solve for elastic GF
ht = md_ISSM.solidearth.lovenumbers.h;
lt = md_ISSM.solidearth.lovenumbers.l;
kt = md_ISSM.solidearth.lovenumbers.k;
[md_GF, vlm_total, vlm_elastic, hlm, accm] = run_greensFunction(md_ISSM, ht, lt, kt, false, lat_gnss, lon_gnss);

% process gnss data
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

% Get VLM from the ISSM solution
[~, ~, ~, ~, delta_B, ~, time_md_pres] = get_GIAoutputs(md_ISSM);

disp('Interpolating present-day model data field onto GPS sites')
interped_matrix = interp_model_to_gnss(md_ISSM, delta_B, lat_gnss, lon_gnss, 'nearest');

temp_interp = false;
for n = 1:length(data_gnss)
    if temp_interp
        % Interpolate onto common times
        interped_md_data_ISSM{n} = interp1(time_md_pres, interped_matrix(n,:), common_time{n}, 'linear', 'extrap');
    else
        % No temporal interpolation
        interped_md_data_ISSM{n} = interped_matrix(n,:); % use this if you ever want to use the original model data
    end

    % Remove mean (no need to detrend linearly because no GIA signal is included)
    disp('Removing mean from the interpolated data')
    interped_md_data_ISSM{n} = detrend(interped_md_data_ISSM{n}, 0);
    % Convert to mm, consistent with the GNSS data
    interped_md_data_ISSM{n} = interped_md_data_ISSM{n} * 1e3;
end

% Get VLM from the GF solution
for n = 1:length(data_gnss)
    if temp_interp
        interped_md_data_GF{n} = interp1(time_md_pres, vlm_total(n,:), common_time{n}, 'linear', 'extrap');
    else
        interped_md_data_GF{n} = vlm_total(n,:);
    end

    disp('Removing mean from the interpolated data')
    interped_md_data_GF{n} = detrend(interped_md_data_GF{n}, 0);

    % Convert to mm
    interped_md_data_GF{n} = interped_md_data_GF{n} * 1e3;
end

%% ==========================================================
%  Multi-panel comparison plot (5 x 4 per figure, all stations)
% ==========================================================
nStations = length(data_gnss);
nRows = 5; 
nCols = 4;
plotsPerFigure = nRows * nCols;

for figIdx = 1:ceil(nStations / plotsPerFigure)
    figure('Name', sprintf('VLM Comparison GF vs. ISSM (Page %d)', figIdx), ...
           'Color','w', 'Position',[100 100 1400 900]);
    tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

    % Determine which stations to plot in this figure
    startIdx = (figIdx-1)*plotsPerFigure + 1;
    endIdx   = min(figIdx*plotsPerFigure, nStations);

    for n = startIdx:endIdx
        nexttile;
        hold on;

        % Plot model results
        if temp_interp
            plot(common_time{n}, interped_md_data_GF{n}, 'b', 'LineWidth', 1.1, 'DisplayName', 'GF');
            plot(common_time{n}, interped_md_data_ISSM{n}, 'r', 'LineWidth', 1.1, 'DisplayName', 'ISSM');
            if ~isempty(data_gnss{n})
                plot(common_time{n}, data_gnss{n}, 'k', 'LineWidth', 1, 'DisplayName', 'GNSS');
            end
            xlabel('Time (year)');
        else
            plot(interped_md_data_GF{n}, 'b', 'LineWidth', 1.1, 'DisplayName', 'GF');
            plot(interped_md_data_ISSM{n}, 'r', 'LineWidth', 1.1, 'DisplayName', 'ISSM');
            xlabel('Timestep');
        end

        % Labels, formatting
        ylabel('meters');
        title(sprintf('Station %s', stn_id{n}), 'FontWeight','bold', 'FontSize', 10);
        grid on;
        set(gca, 'FontSize', 9);
        axis tight;

        % Show legend only in first subplot of each figure
        if n == startIdx
            legend('show', 'Location', 'best', 'FontSize', 8);
        end
    end

    % Title for the page
    sgtitle('Detrended VLM Comparison GF vs. ISSM', 'FontWeight','bold', 'FontSize', 14);
end

% %% ==========================================================
%  Quantitative comparison between ISSM and GF results
%  ==========================================================
disp('Computing RMS difference and correlation between ISSM and GF VLM solutions...');

% Preallocate
nStations = length(data_gnss);
RMS_err = zeros(nStations,1);
corr_coeff = zeros(nStations,1);

for n = 1:nStations
    % Ensure same length vectors (some safety trimming)
    v1 = interped_md_data_ISSM{n}(:);
    v2 = interped_md_data_GF{n}(:);
    len = min(length(v1), length(v2));
    v1 = v1(1:len);
    v2 = v2(1:len);

    % Compute RMS and correlation
    RMS_err(n) = sqrt(mean((v1 - v2).^2, 'omitnan'));
    corr_coeff(n) = corr(v1, v2, 'Rows', 'complete');
end

% Print summary table
fprintf('\n=== GF vs ISSM Vertical Land Motion Agreement ===\n');
fprintf(' %-8s | %-12s | %-10s\n', 'Station', 'RMS diff (mm)', 'Corr coeff');
fprintf('---------------------------------------------\n');
for n = 1:nStations
    fprintf(' %-8s | %12.3f | %10.3f\n', stn_id{n}, RMS_err(n), corr_coeff(n));
end
fprintf('---------------------------------------------\n');
fprintf(' Mean RMS = %.3f mm | Mean Corr = %.3f\n', mean(RMS_err), mean(corr_coeff));
fprintf('=============================================\n\n');

%% ----------------------------------------------------------
% Visualization of agreement metrics
% ----------------------------------------------------------
figure('Name','ISSM vs GF Agreement Summary','Color','w','Position',[200 200 1000 600]);

subplot(2,1,1)
bar(RMS_err)
ylabel('RMS difference (mm)')
title('Greens Function vs ISSM VLM RMS Difference')
set(gca,'XTick',1:nStations,'XTickLabel',stn_id)
grid on

subplot(2,1,2)
bar(corr_coeff)
ylabel('Correlation coefficient')
ylim([0 1])
title('Greens Function vs ISSM VLM Correlation')
set(gca,'XTick',1:nStations,'XTickLabel',stn_id)
xlabel('GNSS Station')
grid on

sgtitle('Benchmarking of Elastic GF against ISSM solutions','FontWeight','bold');

