function correct_ice_loading = adjustIceLoadingForSphere(lat, ice_loading_ellipsoid)
    % Define WGS84 Ellipsoid parameters
    a = 6378137.0;  % Equatorial radius (in meters)
    b = 6356752.3;  % Polar radius (in meters)
    e = sqrt(1 - (b^2 / a^2));  % Eccentricity of the ellipsoid

    % Sphere radius (assuming spherical Earth)
    R_sphere = 6371000;  % Average radius of the Earth in meters (can be adjusted as needed)

    % Calculate the ellipsoid radius of curvature for each latitude
    r_ellipsoid = a * (1 - e^2) ./ (1 - e^2 * sind(lat).^2).^(3/2);

    % Compute the area ratio between ellipsoid and sphere
    area_ratio = (r_ellipsoid ./ R_sphere).^2;

    % Correct the ice loading based on the area ratio
    ice_loading_corrected = ice_loading_ellipsoid .* area_ratio;

    % Output the corrected ice loading
    disp('Ice loading corrected for sphere:');
    disp(ice_loading_corrected);
end
