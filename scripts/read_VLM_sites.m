%% read_VLM_sites.m
% Reads 3D and 1D Vertical Land Motion (VLM) data at GNET GNSS sites in Greenland.
% Optionally matches numbered sites to station names using a coordinate list file
% (e.g., 'station_coords.txt') that follows:
%   KELY   66.98741821   -50.94483762
%
% Usage:
%   data = read_VLM_sites('VLM_sites.dat', true);

function data = read_VLM_sites(filename, plot_fig)

    if nargin < 1
        filename = 'VLM_sites.dat';
    end
    if nargin < 2
        plot_fig = true;
    end

    % --- Read VLM data ---
    opts = detectImportOptions(filename, 'FileType', 'text');
    opts.Delimiter = {'\t',' ',','};
    opts.VariableNames = {'Site','Longitude','Latitude',...
                          'VLM3D_mean','VLM3D_sigma',...
                          'VLM1D_mean','VLM1D_sigma'};
    opts.VariableTypes = {'double','double','double',...
                          'double','double','double','double'};
    tbl = readtable(filename, opts);

    % --- Pack into structure ---
    data.site        = tbl.Site;
    data.lon         = tbl.Longitude;
    data.lat         = tbl.Latitude;
    data.VLM3D_mean  = tbl.VLM3D_mean;
    data.VLM3D_sigma = tbl.VLM3D_sigma;
    data.VLM1D_mean  = tbl.VLM1D_mean;
    data.VLM1D_sigma = tbl.VLM1D_sigma;

    % --- Print summary ---
    fprintf('\n✅ Successfully read %d GNET sites from %s\n', height(tbl), filename);
    fprintf('   Longitude range: %.2f° to %.2f°\n', min(data.lon), max(data.lon));
    fprintf('   Latitude range : %.2f° to %.2f°\n', min(data.lat), max(data.lat));
    fprintf('   Mean 1D VLM    : %.3f ± %.3f mm/yr\n', mean(data.VLM1D_mean), std(data.VLM1D_mean));
    fprintf('   Mean 3D VLM    : %.3f ± %.3f mm/yr\n', mean(data.VLM3D_mean), std(data.VLM3D_mean));

    %% --- Match station names (using same logic as preprocess_gnss_data) ---
    run('settings_observation_data.m');
    coordfile = fullfile(fpath_gnss, fname_coord_gnss);
    if isfile(coordfile)
        fprintf('\n🔍 Matching station names using %s ...\n', coordfile);
        coord_lines = readlines(coordfile);
        coord_lines = coord_lines(~startsWith(coord_lines, '%') & strlength(strtrim(coord_lines)) > 0);

        station_names = {};
        station_lat = [];
        station_lon = [];
        for i = 1:length(coord_lines)
            parts = strsplit(strtrim(coord_lines(i)));
            if numel(parts) >= 3
                station_names{end+1,1} = parts{1};
                station_lat(end+1,1) = str2double(parts{2});
                station_lon(end+1,1) = str2double(parts{3});
            end
        end

        nSites = height(tbl);
        data.name = strings(nSites,1);
        data.match_distance_km = nan(nSites,1);

        for i = 1:nSites
            dlat = data.lat(i) - station_lat;
            dlon = (data.lon(i) - station_lon) .* cosd(data.lat(i));
            dist_km = sqrt((111*dlat).^2 + (111*dlon).^2);
            [mindist, idx] = min(dist_km);
            data.name(i) = station_names{idx};
            data.match_distance_km(i) = mindist;
        end

        % --- Print mapping summary ---
        fprintf('\nMatched stations:\n');
        fprintf('----------------------------------------------\n');
        fprintf('%4s  %-6s  %8s  %8s  %8s\n','ID','Name','Lat(°)','Lon(°)','Δdist(km)');
        fprintf('----------------------------------------------\n');
        for i = 1:nSites
            fprintf('%4d  %-6s  %8.3f  %8.3f  %8.2f\n',...
                data.site(i), data.name(i), data.lat(i), data.lon(i), data.match_distance_km(i));
        end
        fprintf('----------------------------------------------\n');
        fprintf('Average name match distance: %.2f km\n\n', mean(data.match_distance_km));

    else
        fprintf('\n⚠️ No station_coords.txt file found. Skipping name matching.\n\n');
        data.name = repmat("", height(tbl), 1);
        data.match_distance_km = NaN(height(tbl), 1);
    end

    if plot_fig
        % === Create red-white-blue diverging colormap ===
        cmap = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1);   % blue→white
                ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];  % white→red

        % Helper inline for symmetric color axis
        set_caxis = @(vals) caxis([-max(abs(vals(:))) max(abs(vals(:)))]);

        % --- Main figure: 1D and 3D panels side by side ---

        % --- 1D (vertical) VLM ---
        figure('Color','w','Position',[400 200 900 550]);
        scatter(data.lon,data.lat,80,data.VLM1D_mean,'filled');
        hold on;
        if any(data.name~="")
            text(data.lon,data.lat,data.name,'FontSize',8,...
                'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
        colormap(gca,cmap);
        set_caxis(data.VLM1D_mean);
        cb = colorbar;
        cb.Label.String = 'VLM (mm/yr)';
        cb.Label.FontSize = 12;
        cb.Position = [0.92 0.25 0.01 0.55]; 
        cb.Position = cb.Position + [0.03 0 0 0]; 
        xlabel('Longitude (°)');
        ylabel('Latitude (°)');
        title('VLM (1D mean)','FontSize',13,'FontWeight','bold');
        axis equal tight;
        set(gca,'FontSize',12,'LineWidth',1,'Position',get(gca,'Position')+[0 -0.02 0 0.03]); % nudge down slightly

        % --- 3D total VLM ---
        figure('Color','w','Position',[400 200 900 550]);
        scatter(data.lon,data.lat,80,data.VLM3D_mean,'filled');
        hold on;
        if any(data.name~="")
            text(data.lon,data.lat,data.name,'FontSize',8,...
                'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
        colormap(gca,cmap);
        set_caxis(data.VLM3D_mean);
        cb = colorbar;
        cb.Label.String = 'VLM (mm/yr)';
        cb.Label.FontSize = 12;
        cb.Position = [0.92 0.25 0.01 0.55]; 
        cb.Position = cb.Position + [0.03 0 0 0]; 
        xlabel('Longitude (°)');
        ylabel('Latitude (°)');
        title('VLM (3D mean)','FontSize',13,'FontWeight','bold');
        axis equal tight;
        set(gca,'FontSize',12,'LineWidth',1,'Position',get(gca,'Position')+[0 -0.02 0 0.03]);

        % --- Separate figure: 3D − 1D difference ---
        figure('Color','w','Position',[400 200 900 550]);
        diff_vals = data.VLM3D_mean - data.VLM1D_mean;
        scatter(data.lon,data.lat,80,diff_vals,'filled');
        if any(data.name~="")
            text(data.lon,data.lat,data.name,'FontSize',8,...
                'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
        colormap(gca,cmap);
        set_caxis(diff_vals);
        cb = colorbar;
        cb.Label.String = '3D − 1D difference (mm/yr)';
        cb.Label.FontSize = 12;
        cb.Position = [0.92 0.25 0.01 0.55]; 
        cb.Position = cb.Position + [0.03 0 0 0]; 
        xlabel('Longitude (°)');
        ylabel('Latitude (°)');
        title('Difference: 3D − 1D VLM [mm/yr]','FontSize',13,'FontWeight','bold');
        axis equal tight;
        set(gca,'FontSize',12,'LineWidth',1,'Position',get(gca,'Position')+[0 -0.02 0 0.03]);
    end
