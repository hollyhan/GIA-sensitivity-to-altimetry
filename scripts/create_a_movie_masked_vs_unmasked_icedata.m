%% ==============================================================
%  Movie of masked vs. unmasked dh/dt around a GNSS site (e.g., KAGA)
%  Includes earliest & latest mask contours.
% ==============================================================

% --- load data ---
load('dhdt_masked.mat');

% --- Prescribe data to plot ---
data_masked_to_plot = dhdt_masked{1};
data_unmasked_to_plot = dhdt_annual_1;
years_masked = dhdt_masked_years{1};
years_unmasked = years_altimetry_1;
lat = lat_sphere_1; % make sure lat long values are computed from Step 1 of main_script.m
long = long_sphere_1;
dataset_name = 'measureItsLive-GEMB';
md_masked_to_plot = md1_regional_with_mask; % make sure to load this
md_unmasked_to_plot = md1_regional_without_mask;

% --- Color map setup ---
np = 64;
blue_to_white = [linspace(0,1,np/2)', linspace(0,1,np/2)', ones(np/2,1)];
white_to_red  = [ones(np/2,1), linspace(1,0,np/2)', linspace(1,0,np/2)'];
custom_cmap   = [blue_to_white; white_to_red];

% --- Derive mask contours for first & last slices ---
mask_first = dhdt_masked_to_plot(:,:,1);
mask_last  = dhdt_masked_to_plot(:,:,end);
mask_contour_first = ~isnan(mask_first) & mask_first ~= 0;
mask_contour_last  = ~isnan(mask_last)  & mask_last  ~= 0;

% --- Initialize video writer ---
idx_stn = 43;%21; %KUAQ; 43 % KAGA
ylim0 = 69.1;%68.4; %KUAQ %69.1; %KAGA
ylim1 = 69.3;%68.7; %69.3; %KAGA
xlim1 = -49.5;%-33.8; %KUAQ -49.5%KAGA
xlim0 = -50.5; %-32.6; % -50.5;%KAGA
latlim = [ylim0 ylim1];
lonlim = [xlim0 xlim1];
[xmin, ymin] = ll2xy(latlim(1), lonlim(1), +1);
[xmax, ymax] = ll2xy(latlim(2), lonlim(2), +1);

% Convert GNSS site from lat/lon to projection coordinates (EPSG:3413)
[x_gnss, y_gnss] = ll2xy(lat_gnss(idx_stn), lon_gnss(idx_stn), +1);  % +1 for north polar stereographic (Greenland)
[x_gnss_whole, y_gnss_whole] = ll2xy(lat_gnss, lon_gnss, +1);  % +1 for north polar stereographic (Greenland)

vname = sprintf('dhdt_mask_vs_unmask_%s_%s_altimetry_grid.mp4', stn_id{idx_stn}, dataset_name);
v = VideoWriter(vname, 'MPEG-4');
v.FrameRate = 2;      % frames per second
v.Quality   = 100;    % best quality
open(v);

%% --- Loop through time slices ---
for kk = 1:size(data_masked_to_plot,3)

    % --- Figure setup ---
    figure(1); clf;
    set(gcf, 'Position', [100 100 1500 500]);
    set(gcf, 'Renderer', 'zbuffer');   % stable, non-GPU renderer
    set(gcf, 'Color', 'w');

    % ---------- Panel 1: Unmasked ----------
    subplot(1,3,1)
    B = data_unmasked_to_plot(:,:,kk);
    B(B == 0) = NaN;
    pcolor(long, lat, B); shading flat;
    hold on
    %contour(long, lat, mask_contour_first, [1 1], 'g', 'LineWidth', 1);
    %contour(long, lat, mask_contour_last,  [1 1], 'y', 'LineWidth', 1);
    plot(lon_gnss(idx_stn), lat_gnss(idx_stn), 'k.', 'MarkerSize', 15);
    hold off
    title(sprintf('Unmasked (%d-%d)', years_unmasked(kk), years_unmasked(kk+1)));
    set(gca, 'FontSize', 12);
    colormap(flip(custom_cmap)); caxis([-20 20]);
    colorbar;
    ylim([ylim0 ylim1]);
    xlim([xlim0 xlim1]);

    % ---------- Panel 2: Masked ----------
    subplot(1,3,2)
    A = data_masked_to_plot(:,:,kk);
    A(A == 0) = NaN;
    pcolor(long, lat, A); shading flat;
    hold on
    %contour(long, lat, mask_contour_first, [1 1], 'g', 'LineWidth', 1);
    %contour(long, lat, mask_contour_last,  [1 1], 'y', 'LineWidth', 1);
    plot(lon_gnss(idx_stn), lat_gnss(idx_stn), 'k.', 'MarkerSize', 15);
    hold off
    title(sprintf('Masked (%d-%d)', years_masked(kk), years_masked(kk+1)));
    set(gca, 'FontSize', 12);
    colormap(flip(custom_cmap)); caxis([-20 20]);
    colorbar;
    ylim([ylim0 ylim1]);
    xlim([xlim0 xlim1]);

    % ---------- Panel 3: Difference ----------
    subplot(1,3,3)
    C = data_unmasked_to_plot(:,:,kk) - data_masked_to_plot(:,:,kk);
    C(C == 0) = NaN;
    pcolor(long, lat, C); shading flat;
    hold on
    %contour(long, lat, mask_contour_first, [1 1], 'g', 'LineWidth', 1);
    %contour(long, lat, mask_contour_last,  [1 1], 'y', 'LineWidth', 1);
    plot(lon_gnss(idx_stn), lat_gnss(idx_stn), 'k.', 'MarkerSize', 15);
    hold off
    title('Diff in dhdt: unmasked minus masked');
    set(gca, 'FontSize', 12);
    colormap(flip(custom_cmap)); caxis([-20 20]);
    colorbar;
    ylim([ylim0 ylim1]);
    xlim([xlim0 xlim1]);

    % --- Export frame (safe, no OpenGL glitch) ---
    drawnow;
    exportgraphics(gcf, 'temp_frame.png', 'Resolution', 200);
    frame = imread('temp_frame.png');
    writeVideo(v, frame);
end

%% --- Wrap up ---
close(v);
if exist('temp_frame.png','file'); delete('temp_frame.png'); end
disp(['Movie saved as: ', vname]);
%===========================================================================================

% 2. === Draw video on the ISSM mesh ====
vname = sprintf('dhdt_mask_vs_unmask_%s_%s_ISSM_mesh.mp4', stn_id{idx_stn}, dataset_name);
v = VideoWriter(vname, 'MPEG-4');
v.FrameRate = 2;      % frames per second
v.Quality   = 100;    % best quality
open(v);

%% --- Loop through time slices ---
for kk = 1:size(data_masked_to_plot,3)

    % Compute differences between masked vs unmasked ISSM outputs
    diff_h_refined_masked   = md_masked_to_plot.masstransport.spcthickness(1:end-1,kk+1) ...
                            - md_masked_to_plot.masstransport.spcthickness(1:end-1,kk);
    diff_h_refined_unmasked = md_unmasked_to_plot.masstransport.spcthickness(1:end-1,kk+1) ...
                            - md_unmasked_to_plot.masstransport.spcthickness(1:end-1,kk);
    diff_h_refined_diff     = diff_h_refined_unmasked - diff_h_refined_masked;

    % === Use ISSM’s built-in 3-panel plotting ===
    plotmodel(md_unmasked_to_plot, ...
        'data', diff_h_refined_unmasked, ...
        'title', sprintf('Unmasked (%d - %d)',  md_unmasked_to_plot.masstransport.spcthickness(end,kk), md_unmasked_to_plot.masstransport.spcthickness(end,kk+1)), ...
        'data', diff_h_refined_masked, ...
        'title', sprintf('Masked (%d - %d)',  md_masked_to_plot.masstransport.spcthickness(end,kk), md_masked_to_plot.masstransport.spcthickness(end,kk+1)), ...
        'data', diff_h_refined_diff, ...
        'title', 'Unmasked - Masked', ...
        'caxis#all', [-20 20], ...
        'xlim#all', [xmin xmax], ...
        'ylim#all', [ymin ymax], ...
        'colormap#all', flip(custom_cmap), ...
        'colorbar#all', 'on', ...
        'colorbartitle#all', '(m/yr)', ...
        'axis#all', 'equal off', ...
        'fontsize#all', 12);

    % --- Overlay GNSS points (on all panels) ---
    hold on
    plot(x_gnss, y_gnss, 'k.', 'MarkerSize', 15);
    for i = 1:length(stn_id)
        text(x_gnss_whole(i), y_gnss_whole(i), stn_id{i}, ...
             'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r', ...
             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    end

    % --- Write each frame to the video ---
    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% --- Close video writer ---
close(v);


%===========================================================================================
% === If just want snapshots on the ISSM mesh ====
% Differences between masked vs unmasked on the ISSM mesh
diff_h_refined_masked = md_masked_to_plot.masstransport.spcthickness(1:end-1,kk+1)-md_masked_to_plot.masstransport.spcthickness(1:end-1,kk);
diff_h_refined_unmasked =  md_unmasked_to_plot.masstransport.spcthickness(1:end-1,kk+1)-md_unmasked_to_plot.masstransport.spcthickness(1:end-1,kk);

% --- Figure setup ---
figure(1); clf;
set(gcf, 'Position', [100 100 1500 500]);
set(gcf, 'Renderer', 'zbuffer');   % stable, non-GPU renderer
set(gcf, 'Color', 'w');

% ---------- Panel 1: Unmasked ----------
subplot(1,3,1)

plotmodel(md_unmasked_to_plot,'data', diff_h_refined_unmasked,'caxis',([-20 20]), 'xlim',[xmin xmax],'ylim',[ymin ymax])
hold on
title(sprintf('Unmasked (%d-%d)', md_unmasked_to_plot.masstransport.spcthickness(end,kk), md_unmasked_to_plot.masstransport.spcthickness(end,kk+1)))
colormap(flip(custom_cmap));  
plot(x_gnss, y_gnss, 'k.', 'MarkerSize', 15);

for i = 1:length(stn_id)
    text(x_gnss_whole(i), y_gnss_whole(i), stn_id{i}, ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', 'r', ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center');
end

% ---------- Panel 2: Masked ----------
subplot(1,3,2)
plotmodel(md_masked_to_plot,'data', diff_h_refined_masked,'caxis',([-20 20]), 'xlim',[xmin xmax],'ylim',[ymin ymax])
hold on
title(sprintf('Masked (%d-%d)',md_masked_to_plot.masstransport.spcthickness(end,kk),md_masked_to_plot.masstransport.spcthickness(end,kk+1)))
colormap(flip(custom_cmap));
plot(x_gnss, y_gnss, 'k.', 'MarkerSize', 15);

for i = 1:length(stn_id)
    text(x_gnss_whole(i), y_gnss_whole(i), stn_id{i}, ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', 'r', ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center');
end

% ---------- Panel 3: Difference ----------
subplot(1,3,3)
plotmodel(md_masked_to_plot,'data', diff_h_refined_unmasked-diff_h_refined_masked,'caxis',([-20 20]), 'xlim',[xmin xmax],'ylim',[ymin ymax],'edgecolor','black')
hold on
title(sprintf('Unmasked minus Masked (%d-%d)'))
colormap(flip(custom_cmap));  
plot(x_gnss_whole, y_gnss_whole, 'm.', 'MarkerSize', 15);
plot(x_gnss, y_gnss, 'b.', 'MarkerSize', 15);

for i = 1:length(stn_id)
    text(x_gnss_whole(i), y_gnss_whole(i), stn_id{i}, ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', 'r', ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center');
end
%===========================================================================================

