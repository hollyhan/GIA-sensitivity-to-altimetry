clear; clc; close all;

% Process altimetry datasets without firn correction
disp('=== Processing altimetry data ===');
[h_annual_1, dhdt_annual_1, dhdt_monthly_1, years_altimetry_1, lat_sphere_1, long_sphere_1] = preprocess_ice_altimetry('measureItsLive');
[h_annual_2, dhdt_annual_2, dhdt_monthly_2, years_altimetry_2, lat_sphere_2, long_sphere_2] = preprocess_ice_altimetry('DTU2016');
[h_annual_3, dhdt_annual_3, dhdt_monthly_3, years_altimetry_3, lat_sphere_3, long_sphere_3] = preprocess_ice_altimetry('DTU2025');
disp('====================================');

% Process firn model datasets
disp('=== Processing firn model data ===');
[dfac_annual_1, dfac_monthly_1, years_firn_1, lat_sphere_firn_1, long_sphere_firn_1] = preprocess_firn_model('GEMB');
[dfac_annual_2, dfac_monthly_2, years_firn_2, lat_sphere_firn_2, long_sphere_firn_2] = preprocess_firn_model('GSFC');
[dfac_annual_3, dfac_monthly_3, years_firn_3, lat_sphere_firn_3, long_sphere_firn_3] = preprocess_firn_model('RACMO2.3p2');
disp('====================================');

% Process glacier mask datasets (returns native high-resolution mask on spherical geographic coordinates)
disp('=== Processing glacier mask data ===');
[ice_masks, mask_years, lat_mask, lon_mask, x_mask, y_mask] = preprocess_glacier_mask([1996:1997]);
disp('====================================');

% Resample glacier mask to each altimetry dataset's EXACT native grid
disp('=== Resampling glacier mask to each altimetry dataset ===');
ice_masks_measureItsLive = resample_mask_to_target_grid(ice_masks, lat_mask, lon_mask, lat_sphere_1, long_sphere_1);
ice_masks_DTU2016 = resample_mask_to_target_grid(ice_masks, lat_mask, lon_mask, lat_sphere_2, long_sphere_2);
ice_masks_DTU2025 = resample_mask_to_target_grid(ice_masks, lat_mask, lon_mask, lat_sphere_3, long_sphere_3);
disp('====================================');

