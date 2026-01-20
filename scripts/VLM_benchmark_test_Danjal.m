% Read in the two ice altimetry data files
% Data 1.
disp("Loading ice elevation data used in Berg et al., 2024")

np = 64; % Number of colors
blue_to_white = [linspace(0,1,np/2)', linspace(0,1,np/2)', ones(np/2,1)];
white_to_red = [ones(np/2,1), linspace(1,0,np/2)', linspace(1,0,np/2)'];
custom_cmap = [blue_to_white; white_to_red];

% For the first two datasets (Berg et al., 2024)
% To convert to ice thickness change from water equivalent height 
rhoo = 1000.0; % density of water (kg/m^3)
rhoi =  917.0;  % density of ice (kg/m^3)
water2ice = rhoo / rhoi;
kg_per_Gt = 1e12;
pixel_area = 1000 * 1000;  % 1 km x 1 km in m^2
% Create time vector if needed
years = 2011:2019;   

% --- Read altimetry file ---
%Dataset 1
fname = '/Users/kyhan/Downloads/dhdt_2011-2020_GrIS_1km_AUG2021.txt';

% Read everything as numeric
data = readmatrix(fname);

% Extract latitude and longitude
lat = data(:,1);
lon = data(:,2);

% dh and error come in alternating pairs:
% [dh_2011_2012, err_2011_2012, dh_2012_2013, err_2012_2013, ...]
vals = data(:,3:end);

% Number of years
nyears = size(vals,2) / 2;

% Preallocate dh and sigma
dh     = zeros(size(lat,1), nyears);
sigma  = zeros(size(lat,1), nyears);

% Fill
dh(:, :)    = vals(:, 1:2:end);
sigma(:, :) = vals(:, 2:2:end);


dh(:,:) = dh(:,:) .* water2ice;
sigma(:,:) = sigma(:,:) .* water2ice;
   % Centered on the 15 Apr – 15 Apr interval
% Create projection objects
wgs84 = geocrs(4326);       % WGS84 lat/lon
greenland_ps = projcrs(3413); % NSIDC Sea Ice Polar Stereographic North

% Convert
[x, y] = projfwd(greenland_ps, lat, lon);

mass_change_kg_per_pixel_1 = rhoi * dh * pixel_area;   % [kg/yr]

% Sum over all pixels
mass_change_total_kg_per_year_1 = nansum(mass_change_kg_per_pixel_1, 1);

% Convert to Gt/yr
mass_change_Gt_per_year_1 = mass_change_total_kg_per_year_1 / kg_per_Gt;

figure('Color','w');
scatter(x, y, 6, dh(:,1), 'filled');
axis equal tight;
colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('dh/dt (m/yr) 2011–2012 (native 1 km grid, EPSG:3413)');
colormap(flip(custom_cmap));   % or your preferred colormap
caxis([-1 1])
plot(years, mass_change_Gt_per_year,'ro') 
hold on
plot(years, mass_change_Gt_per_year) 
ylim([-1000 1000])

 =============================================================================================
% Data 2. 
fname2 = '/Users/kyhan/Downloads/dhdt_2011-2020_GrIS_1km_final_Berg_et_al.parquet';
T = parquetread(fname2);

years = 2011:2019;

for i = 1:length(years)
    field = sprintf('x%d', years(i));
    T.(field) = T.(field) * water2ice;
end

ps = projcrs(3413);
[T.x_ps, T.y_ps] = projfwd(ps, T.lat, T.lon);

scatter(T.x_ps, T.y_ps, 6, T.x2011, 'filled');
colorbar
axis equal tight
colormap(flip(custom_cmap));   % or your preferred colormap
caxis([-1 1])

% Calculate mass change GT/year
mass_change_Gt_per_year_2 = zeros(1, length(years));
for i = 1:length(years)
    field = sprintf('x%d', years(i));   % e.g., 'x2011'
    
    dh_ice = T.(field);                 % m ice / yr
    
    % Convert to mass per pixel (kg)
    mass_change_pixel_kg_2 = rhoi * dh_ice * pixel_area;
    
    % Total mass change (kg)
    mass_total_kg_2 = nansum(mass_change_pixel_kg_2);
    
    % Convert to Gt
    mass_change_Gt_per_year_2(i) = mass_total_kg_2 / kg_per_Gt;
end


% Based on the above, we can conclude that the two files yield the very similar amount of mass loss.

% Data 3. =============================================================================================
disp("Loading ice elevation data from DTU2025")
filename = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/altimetry/khan_et_al_2025/Greenland_dhdt_icevol_1kmgrid_DB.nc';
fname_firn = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/altimetry/khan_et_al_2025/Greenland_dhdt_firn_1kmgrid_DB.nc';
fname_elas = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/altimetry/khan_et_al_2025/Greenland_dhdt_elas_1kmgrid_DB.nc';
fname_GIA = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/altimetry/khan_et_al_2025/Greenland_GIA_1kmgrid_DB.nc';

% Read the variables
X = ncread(filename, 'x');       % East coordinate (m)
Y = ncread(filename, 'y');       % North coordinate (m)
Time = ncread(filename, 'time'); % Time in days since 2003-01-01

% Create a 3D geographic coordinate for the analysis for regional mass conservation the comes later in the script
[Xgrid, Ygrid] = meshgrid(X, Y); %for the DTU2025 data
[lat_dtu, lon_dtu] = projinv(projcrs(3413), Xgrid, Ygrid);  % WGS84 lat/lon
x_dtu2025 = 6371012 .* cosd(lat_dtu) .* cosd(lon_dtu);
y_dtu2025 = 6371012 .* cosd(lat_dtu) .* sind(lon_dtu);
z_dtu2025 = 6371012 .* sind(lat_dtu);

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

% Calculate the total monthly change in mass (m/month)
dhdt_monthly = dhdt_vol - dhdt_firn + dhdt_elas + dhdt_GIA;

% Convert fill values to NaN (redundant but to be safe)
dhdt_monthly(dhdt_monthly == -9999) = NaN;

% Convert time from days since 2003-01-01 to decimal years
reference_date = datetime(2003, 1, 1);
Time_datetime = reference_date + days(Time);
Time = year(Time_datetime) + (day(Time_datetime, 'dayofyear') - 1) / 365.25;

% Filter data to end of 2022 (exclude incomplete 2023 data)
time_mask = (Time >= 2011) & (Time <= 2020); % Include only data between 2011-2020
Time = Time(time_mask);
dhdt_monthly = dhdt_monthly(:, :, time_mask);

% Define years vector for output from actual Time data
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

dhdt_total = sum(dhdt_annual, 3, 'omitnan');

% calcualte mass change in Gt per year 
dx = diff(X(1:2));
dy = diff(Y(1:2));
pixel_area = dx * dy;

mass_change_Gt_per_year_3 = zeros(1, length(years));
mass_change_kg_per_pixel_3  = zeros(length(X), length(Y), length(years)); % this is for the later step of calculating regional mass change
for i = 1:length(years)
    dh_ice = dhdt_annual(:,:,i);              % m ice / yr
    
    % Convert to mass per pixel (kg)
    mass_change_kg_per_pixel_per_year_3 = rhoi .* dh_ice .* pixel_area;
    mass_change_kg_per_pixel_3(:,:,i) = mass_change_kg_per_pixel_per_year_3;
    
    % Total mass change (kg)    
    mass_total_kg_3 = squeeze(sum(mass_change_kg_per_pixel_per_year_3, [1, 2], 'omitnan')); %

    % Convert to Gt
    mass_change_Gt_per_year_3(i) = mass_total_kg_3 / kg_per_Gt;
end

% plot ==========
figure('Color','w', 'Position', [100 100 1500 500]);
subplot(1,2,1);
plot(years, -mass_change_Gt_per_year_1, 'ko', 'DisplayName','Berg et al., text format');
hold on
plot(years, -mass_change_Gt_per_year_2, 'r', 'DisplayName','Berg et al., Parquet format');
plot(years, -mass_change_Gt_per_year_3, 'bo-', 'DisplayName', 'DTU2025')
xlabel('Year');
ylabel('Mass Loss (Gt/yr)');
title('Mass Loss between 2011-2020 from different DTU sources (on native grid)');
grid on
legend
set(gca, 'FontSize', 16); 

subplot(1,2,2);
plot(years, -(mass_change_Gt_per_year_3 - mass_change_Gt_per_year_2), 'ko','MarkerFaceColor', 'k');
hold on
xlabel('Year');
ylabel('Mass change (Gt/yr)');
title('Differences in Mass loss (DTU2025 minus Berg et al., 2024 data)');
grid on
yline(0,'k-','LineWidth',1);
set(gca, 'FontSize', 16); 

% Check mass variations between 2011-2022 for the area around GNSS stations (e.g., ASKY and HEL2 stations) from the two datasets
% on the native grid

% Read in the GNSS station information
run('settings_observation_data.m');
[lat_gnss, lon_gnss, data_gnss, err_gnss, time_gnss, R2_gnss] = preprocess_gnss_data(stn_id, fpath_gnss, fname_coord_gnss, n_degree);

% calculate ice mass variations within a certain distance (e.g., 10km) to gnss stations
% I have, (1) lat and long location of each gnss station, (2) x and y coordinates (and lat and long) of ice altimetry data points,
% (3) mass change per year in Gt per altimerty pixel. So what needs to be done is 
% for each lat and long location of GNSS data, select pixels that are within a chosen distance 'dist', and sum the mass cange per year in for that radius.
% Since the first two datasets are very similar, just use one of them.

distance = 100000 % distance in meters
% Get x and y positions of lat_pt and lon_pt
greenland_ps = projcrs(3413); % NSIDC Sea Ice Polar Stereographic North
% Convert
[x_gnss, y_gnss] = projfwd(greenland_ps, lat_gnss, lon_gnss);
[Xgrid, Ygrid] = meshgrid(X, Y); %for the DTU2025 data
%create arrays for each region around gnss stations
dmass_region_1 = zeros(length(stn_id), length(years));
dmass_region_2 = zeros(length(stn_id), length(years));

for n = 1:length(stn_id)
    % For dataset 1
    x_dist = x_gnss(n) - x; 
    y_dist = y_gnss(n) - y;
    dist_pt = sqrt(x_dist.^2 + y_dist.^2);
    % find pixels that lie within dist
    idx = find(dist_pt(dist_pt <= distance));

    for t = 1:length(years)
        dmass_region_1(n,t) = nansum(mass_change_kg_per_pixel_1(idx,t), 1) / kg_per_Gt;
    end

    % For dataset 3 (DTU2025)
    x_dist = x_gnss(n) - Xgrid;
    y_dist = y_gnss(n) - Ygrid;
    dist_pt = sqrt(x_dist.^2 + y_dist.^2); % This doesnt work because of different sizes of X and Y. How do I get it work?
    
    % find pixels that lie within dist
    idx = (dist_pt <= distance);  % logical mask same size as grid

     % Sanity check: area of masked DTU region
    area_region_native(n) = nnz(idx) * pixel_area;

    % mass change per pixel: 2520 × 1462 × Nyears
    % sum over spatial dimension using mask
    for t = 1:length(years)
        layer = mass_change_kg_per_pixel_3(:,:,t);
        dmass_region_3(n,t) = nansum(layer(idx))/ kg_per_Gt; % mass change within the distance per Gt
    end
end

% Calcualte differences in trend
for n = 1:length(stn_id)
    p1 = polyfit(years, dmass_region_1(n,:), 1);   % slope of DTU2016
    p3 = polyfit(years, dmass_region_3(n,:), 1);   % slope of DTU2025

    mass_trend_1(n) = p1(1); %slope of local mass change (Gt/yr^2) around station n 
    mass_trend_3(n) = p3(1);
    diff_local_mass_trend(n) = mass_trend_3(n) - mass_trend_1(n);
end

%plot
n = 47 ; %ASKY station

figure('Color','w');
hold on;
% Original curves
plot(years, dmass_region_1(n,:), 'b-o', 'LineWidth', 1.5, ...
    'DisplayName','DTU old');
plot(years, dmass_region_3(n,:), 'r-o', 'LineWidth', 1.5, ...
    'DisplayName','DTU2025');

% ---- Compute linear fits ----
p1 = polyfit(years, dmass_region_1(n,:), 1); % [slope, intercept]
p3 = polyfit(years, dmass_region_3(n,:), 1);
% Trend line over same x-range
years_fit = linspace(min(years), max(years), 100);
trend1 = polyval(p1, years_fit);
trend3 = polyval(p3, years_fit);
plot(years_fit, trend1, 'b--', 'LineWidth', 2, ...
    'DisplayName', sprintf('DTU old trend (%.3f Gt/yr^2)', p1(1)));
plot(years_fit, trend3, 'r--', 'LineWidth', 2, ...
    'DisplayName', sprintf('DTU2025 trend (%.3f Gt/yr^2)', p3(1)));
xlabel('Year');
ylabel('Mass change (Gt/yr)');
title('Mass trend (Gt/yr) for every year with fitted trends at stn:', stn_id{n});
grid on;
legend('Location','best');
set(gca,'FontSize',14);

figure('Color','w','Position',[100 100 1700 600]);
x_pos = 1:length(stn_id);
hold on
%bar(x_pos, diff_local_mass_trend);
bar(x_pos, mass_trend_1);
%bar(x_pos, mass_trend_3);
xlabel('Station ID', 'FontSize', 14);
ylabel('Diff in rate of mass loss between 2011-2020 (Gt/yr^2)', 'FontSize', 16);
titlestr = sprintf('Diff in mass loss rate between DTU-Danjal and DTU2025 within %d km of distance from GNSS sites',distance/1000)
title(titlestr)
grid on;
set(gca, 'FontSize', 12);
xticks(x_pos);
xticklabels(stn_id);
xtickangle(45);

diff_mass_model = dmass_region_3 - dmass_region_1; % Gt/yr difference per year 


%%%%%====================================================
% Load the interpolated ice on the ISSM mesh with the DTU dta
md_dim = '3D_global'; 
if strcmp(md_dim,'2D_regional')
    md = loadmodel('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results/mesh/md4_regional_without_mask.mat');
    md_regional_2D = md;
elseif strcmp(md_dim, '3D_regional')
    md_regional_2D = loadmodel('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results/mesh/md4_regional_without_mask.mat');
    md = TransformCoord_2Dto3D(md_regional_2D, 3413);
    md_regional_3D = md;
elseif strcmp(md_dim, '3D_global')
    md = loadmodel('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results/md4_solved_Elastic_with_mask.mat');
    md_global_3D = md;
end

load('/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results/mass_report_with_mask.mat'); % this one is for the masked version
%plot(2011:2019,mass_report_4_with_mask.mass_src_total(9:17),'b')
%hold on
%plot(2011:2019,mass_report_4_with_mask.mass_dst_total(9:17),'bo')

% Calculate how much mass is being lost over Greenland on the ISSM mesh
% Get size of each element to get total mass change within a cell
num_vert = md.mesh.numberofvertices;
num_elem = md.mesh.numberofelements;
index_elem = md.mesh.elements;
x_vert = md.mesh.x;
y_vert = md.mesh.y;
x_elem = mean(x_vert(index_elem), 2);
y_elem = mean(y_vert(index_elem), 2);

% Calculate element areas
if strcmp(md_dim, '2D_regional');
    a_elem = GetAreas(index_elem, x_vert, y_vert); % m²
elseif strcmp(md_dim,'3D_global') || strcmp(md_dim,'3D_regional')
    z_vert = md.mesh.z;
    z_elem = mean(z_vert(index_elem), 2);
    a_elem = GetAreas3DTria(index_elem, x_vert, y_vert, z_vert); % m²
end

% Calulate mass changes for each cell
rho = md.materials.rho_ice;  % kg/m³
dh_vert = diff(md.masstransport.spcthickness(1:end-1,:), 1, 2);
nt = length(md.masstransport.spcthickness(end,:));
% DEBUG: Check the mass changes
fprintf('=== DEBUGGING MASS CHANGES ===\n');
fprintf('dh_vert matrix size: %dx%d\n', size(dh_vert, 1), size(dh_vert, 2));
fprintf('Total NaN values in dh_vert: %d\n', sum(isnan(dh_vert), 'all'));
fprintf('dh_vert range: %.6f to %.6f\n', min(dh_vert(:)), max(dh_vert(:)));
fprintf('================================\n\n');
% Calculate mass changes per element
% Each column represents instantaneous mass change at that time step
dm = zeros(num_elem, nt-1);
dmass_Gt = zeros(nt-1, 1);
cum_dmass_Gt = zeros(nt-1, 1); % cummulative mass loss across the whole cells
cum_dmass_old = 0;
for n = 1:nt-1
    dh_vert_slice = dh_vert(:, n);
    dh_elem = mean(dh_vert_slice(index_elem), 2, 'omitnan');  % in meters, ignore NaN
    dm(:, n) = rho * (dh_elem .* a_elem);  % kg (mass change per element)

    % DEBUG: Check for NaN in this time step
    if n <= 3  % Only debug first few time steps to avoid spam
        nan_count_slice = sum(isnan(dh_vert_slice));
        nan_count_elem = sum(isnan(dh_elem));
        nan_count_dm = sum(isnan(dm(:, n)));
        fprintf('Time step %d: dh_vert NaN=%d, dh_elem NaN=%d, dm NaN=%d\n', ...
                n, nan_count_slice, nan_count_elem, nan_count_dm);
    end

    dmass_Gt(n) = nansum(dm(:,n)) / kg_per_Gt; ;
    cum_dmass_Gt(n) = cum_dmass_old + nansum(dm(:,n)) / kg_per_Gt;;
    cum_dmass_old = cum_dmass_Gt(n); % cumulative mass loss in Gt
end

if strcmp(md_dim, '2D_regional')
    dmass_Gt_regional_2D = dmass_Gt;
    cum_dmass_Gt_regional_2D = cum_dmass_Gt;
    dm_regional_2D = dm;
    legstr = 'ISSM 2D regional mesh';
elseif strcmp(md_dim, '3D_regional')
    dmass_Gt_regional_3D = dmass_Gt;
    cum_dmass_Gt_regional_3D = cum_dmass_Gt;
    legstr = 'ISSM 3D regional mesh';
elseif strcmp(md_dim, '3D_global')
    dmass_Gt_global_3D = dmass_Gt;
    cum_dmass_Gt_global_3D = cum_dmass_Gt;
    dm_global_3D = dm;
    legstr = 'ISSM 3D global mesh';
end

%%%%%====================================================
% Calculate the Greenland portion of the 3D global mesh 
% Run this step separately from the above three cases. 
md_world = createWorldMesh(md_regional_2D);
[md_global, transfer] = mergeMeshes(md_regional_3D, md_world);
nWorldV = md_world.mesh.numberofvertices;
nRegV  = md_regional_3D.mesh.numberofvertices;
RegV_global = transfer(nWorldV+1 : nWorldV+nRegV, 2); 
RegElems_global = find(any(ismember(md_global.mesh.elements, RegV_global), 2));
elems = md_global.mesh.elements;
index_elem = md_global.mesh.elements; 
x = md_global.mesh.x;
y = md_global.mesh.y;
z = md_global.mesh.z;
a_elem_global = GetAreas3DTria(elems, x, y, z);   % m²
dh_vert = diff(md_global_3D.masstransport.spcthickness(1:end-1,:), 1, 2); 

dm = zeros(num_elem, nt-1);
dmass_Gt_global_3D_GrIS = zeros(nt-1, 1);
cum_dmass_Gt_global_3D_GrIS = zeros(nt-1, 1); % cummulative mass loss across the whole cells
cum_dmass_old = 0;
for n = 1:nt-1
    dh_vert_slice = dh_vert(:, n);
    dh_elem = mean(dh_vert_slice(index_elem), 2, 'omitnan');  % in meters, ignore NaN
    dm(:, n) = rho * (dh_elem .* a_elem_global);  % kg (mass change per element)

    dmass_Gt(n) = nansum(dm(:,n)) / kg_per_Gt; ;
    cum_dmass_Gt_global_3D_GrIS(n) = cum_dmass_old + nansum(dm(:,n)) / kg_per_Gt;;
    cum_dmass_old = cum_dmass_Gt_global_3D_GrIS(n); % cumulative mass loss in Gt
end
dmass_Gt_global_3D_GrIS = dmass_Gt;


% Now, compare mass changes on the native altimetry data, ISSM-2D mesh, 3D regional mesh, and global 3D mesh
% ISSM 2D to 3D regional mesh
% need to have consecutively run the 2D_regional and 3D_regional cases

% After running the 2D_regional, 3D_regional, 3D_global cases
% compare only during the period between 2011-2020
figure('Color', 'w')
%plot(years, -mass_change_Gt_per_year_1, 'k-*', 'DisplayName','DTU2022 on native grid (text format data)');
%hold on
plot(years, -mass_change_Gt_per_year_3, 'r-*', 'DisplayName', 'DTU2025 on native data grid')
hold on
plot(years, -dmass_Gt_regional_2D(9:17), 'b.--', 'DisplayName', 'on ISSM 2D regional mesh')
hold on
plot(years, -dmass_Gt_regional_3D(9:17), 'r-o', 'DisplayName', 'on ISSM 3D regional mesh') 
plot(years, -dmass_Gt_global_3D(9:17), 'k-o', 'DisplayName', 'on ISSM 3D global mesh') 
plot(years, -dmass_Gt_global_3D_GrIS(9:17), 'y-*', 'DisplayName','on ISSM 3D global but GrIS portion only')
% optionally add the 3D global GrIS portion here
legend('Location', 'best')
title('Mass Loss in DTU data between 2011-2020')
xlabel('Years')
ylabel('Mass Loss (Gt/yr)')
grid on

%%%%%====================================================
% Once the Greenland-wide interpolation is confirmed, now check regional mass conservation. 
% Just use the DTU2025 on the native grid and the interpolated version on the 3D ISSM global mesh
% This requires having run the section above that calculates regional mass changes on the DTU native grid

% Convert DTU2025 native data onto 3D_global

distance = 100000; % in meters
%create arrays for each region around gnss stations
nt_global_3D = md_global_3D.masstransport.spcthickness(end,:);
dmass_region_on_ISSM3Dglobal_mesh = zeros(length(stn_id), nt-1);
index_elem = md_global_3D.mesh.elements; 
[xm, ym, zm] = deal(md_global_3D.mesh.x, ...
                    md_global_3D.mesh.y, ...
                    md_global_3D.mesh.z);

lon_mesh = atan2d(ym, xm);
hyp = sqrt(xm.^2 + ym.^2);
lat_mesh = atan2d(zm, hyp);
[Xmesh_3413, Ymesh_3413] = projfwd(projcrs(3413), lat_mesh, lon_mesh);
Xi_elem = mean(Xmesh_3413(index_elem),2);
Yi_elem = mean(Ymesh_3413(index_elem),2);

isGreenlandElem = false(md_global_3D.mesh.numberofelements, 1);
isGreenlandElem(RegElems_global) = true;


R = md_global_3D.solidearth.planetradius;
for n = 1:length(stn_id)
    dist = sqrt( (Xi_elem - x_gnss(n)).^2 + (Yi_elem - y_gnss(n)).^2 );
    elem_mask = dist <= distance;

    local_mask = elem_mask & isGreenlandElem;


    % Compute mass per year inside this region
    for t = 1:nt-1
        dmass_region_on_ISSM3Dglobal_mesh(n,t) = nansum( dm_global_3D(local_mask,t) ) / kg_per_Gt;
        dmass_global_on_ISSM3Dglobal_mesh(t) = nansum(dm_global_3D(:,t)) / kg_per_Gt;
    end
end

%plot
n = 47 ; %ASKY station
figure('Color','w');
hold on;
% Original curves
%plot(years, dmass_region_1(n,:), 'k-o', 'LineWidth', 1.5,'DisplayName','DTU2022 native grid');
plot(years, dmass_region_3(n,:), 'r-o', 'LineWidth', 0.5, 'DisplayName','DTU2025 native grid');
plot(years, dmass_region_on_ISSM3Dglobal_mesh(n,9:17), 'b-*', 'DisplayName','on ISSM 3D global')
legend()
%%% ===============================================
%% 2. Regional mass from ISSM 2D regional mesh

elems2D = md_regional_2D.mesh.elements;     % [nElem2D x 3]
xc2D    = mean(md_regional_2D.mesh.x(elems2D), 2);  % element centers
yc2D    = mean(md_regional_2D.mesh.y(elems2D), 2);

nt = length(md_regional_2D.masstransport.spcthickness(end,:));
dmass_region_ISSM2D = zeros(length(stn_id), nt-1);

for n = 1:length(stn_id)
    dist_elem = hypot(xc2D - x_gnss(n), yc2D - y_gnss(n));
    elem_mask = (dist_elem <= distance);

    fprintf('Station %s: %d 2D elements inside %.0f km\n', ...
            stn_id{n}, nnz(elem_mask), distance/1e3);

    for k = 1:nt-1
        dmass_region_ISSM2D(n,k) = ...
            nansum(dm_regional_2D(elem_mask, k)) / kg_per_Gt;
        test(k) = nansum(dm_regional_2D(:,k))/kg_per_Gt;
    end
end

%%%%%====================================================

% Check for ASKY site
% Back-calculate mean dh using mass change
area_native = area_region_native(47);
for k = 1:nYears
    % DTU: we know regional mass in Gt/yr already
    M_native_Gt = dmass_region_DTU_native(iASKY,k);   % Gt/yr (negative for loss) same as  dmass_region_3
    M_native_kg = M_native_Gt * kg_per_Gt;            % kg/yr
    % Convert to mean dh in m/yr over that disc
    mean_dh_native(k) = M_native_kg / (rho_ice * area_native);
end

area_issm   = nansum(a_elem_2D(elem_mask));
for k = 1:nt-1 %nt = 21
    % ISSM: sum dm over masked elements
    M_issm_Gt =  dmass_region_ISSM2D(iASKY, k);
    M_issm_kg = M_issm_Gt * kg_per_Gt;   % kg/yr
    mean_dh_issm(k) = M_issm_kg / (rho_ice * area_issm);
end

% Direct-calculate mean dh using the dh field
for t = 1:nYears
    % --- DTU native, region mask already defined as mask_native (Ny×Nx) ---
    layer_native = dhdt_annual(:,:,t);     % m/yr on native grid
    dh_native_mean(t) = mean(layer_native(idx),'omitnan');  % m/yr

    % --- ISSM 2D, same 100 km circle using elem_mask_2D (Nelements×1) ---
    dh_issm_elem = dh_elem_2D(:,t);        % m/yr on ISSM elements
    dh_issm_mean(t) = mean(dh_issm_elem(elem_mask_2D),'omitnan');  % m/yr
end

%% Plot comparison
figure('Color','w');
subplot(2,1,1); hold on; grid on;
plot(years, dmass_region_DTU_native(iASKY,:), 'ro-', 'DisplayName','DTU native');
plot(years, squeeze(sum(dm2D_2011_2019(elem_mask,:),1))/kg_per_Gt, ...
     'b*-', 'DisplayName','ISSM 2D');
ylabel('Mass change (Gt/yr)');
title(sprintf('%s, R = %.0f km', stn_id{iASKY}, R/1e3));
legend('Location','best');

subplot(2,1,2); hold on; grid on;
plot(years, mean_dh_native*1e3, 'ro-', 'DisplayName','DTU native');
plot(years, mean_dh_issm*1e3,   'b*-', 'DisplayName','ISSM 2D');
ylabel('Mean dh (mm/yr)');
xlabel('Year');
legend('Location','best');


