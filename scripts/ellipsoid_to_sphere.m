function [lat_sphere, long_sphere, dhdt_sphere] = ellipsoid_to_sphere(lat, long, r_sphere, dhdt_ellipsoid, plot_fig)
    % ellipsoid_to_sphere.m
    % Holly Han (created: Oct. 29th, 2024; Last edited: Nov. 5th, 2024).    
    % Converts latitude, longitude, and height from WGS84 ellipsoid to spherical coordinates
    %
    % Input:
    % lat  - Latitude in degrees
    % long - Longitude in degrees
    % dhdt_ellipsoid - Height changes on the ellipsoid (in meters)
    % r_sphere - Mean radius of the spherical Earth (in meters)
    % plot_fig - (Optional) If true, plots the radius of curvature for debugging
    %
    % Output:
    % lat_sphere  - Latitude on the spherical Earth (unchanged from ellipsoid)
    % long_sphere - Longitude on the spherical Earth (unchanged from ellipsoid)
    % dhdt_sphere - Height relative to Earth's center on the spherical Earth (in meters)

    if nargin < 5
        plot_fig = false;  % Default plot flag to false if not provided
    end

    % Convert latitude and longitude to radians
    lat_rad = deg2rad(lat);

    % Define WGS84 ellipsoid parameters
    a = 6378137.0;            % Semi-major axis (equatorial radius) in meters
    f = 1 / 298.257223563;     % Flattening of WGS84 ellipsoid
    b = a * (1 - f);           % Semi-minor axis (polar radius) in meters

    % Calculate the radius of curvature in the prime vertical at each latitude
    N = a ./ sqrt(1 - (2*f - f^2) * sin(lat_rad).^2);  % Radius of curvature in the prime vertical

    if plot_fig
        % For debugging: Plot the radius of curvature if plot_flag is true
        figure;
        plot(lat, r_sphere ./ N);
        xlabel('Latitude (degrees)');
        ylabel('Radius of Curvature (m)');
        title('Radius of Curvature (N) vs Latitude');
    end

    % Correct for ice thickness change based on curvature
    % Create correction factor with proper dimensions
    correction_factor = r_sphere ./ N;

    % Check if we need to transpose correction_factor
    [nlat_data, nlon_data] = size(dhdt_ellipsoid(:,:,1));
    [nlat_corr, nlon_corr] = size(correction_factor);
    
    if nlat_data == nlat_corr && nlon_data == nlon_corr
        % Dimensions already match, no transpose needed
        disp('Dimensions match, no transpose needed');
    elseif nlat_data == nlon_corr && nlon_data == nlat_corr
        % Need to transpose correction_factor
        correction_factor = correction_factor';
        disp('Transposed correction_factor to match data orientation');
    else
        error('Dimension mismatch that cannot be resolved by transposition');
    end
    
    % Handle 3D data by applying the correction to each time slice
    if ndims(dhdt_ellipsoid) == 3
        [nlat, nlon, ntime] = size(dhdt_ellipsoid);
        dhdt_sphere = zeros(size(dhdt_ellipsoid));
        
        for t = 1:ntime
            % Apply correction to each time slice
            dhdt_sphere(:,:,t) = dhdt_ellipsoid(:,:,t) .* correction_factor;
        end
    else
        % Handle 2D data
        dhdt_sphere = dhdt_ellipsoid .* correction_factor;
    end

    % The latitude and longitude on the spherical Earth remain unchanged
    lat_sphere = lat;
    long_sphere = long;

    % Display a message indicating completion
    disp('Conversion from ellipsoid to spherical coordinates completed.');
end
