function interped_matrix = interp_model_to_gnss(md, field_data, lat_gnss, lon_gnss, method)
% Holly Han (created April 20th, 2025; Last edited in April. 30th, 2025).
% interpolate_model_to_gnss Interpolates a time-varying model field onto GPS sites.
%
% Inputs:
%   md            - Model object from which mesh information is extracted
%   field_data    - [NxT] matrix of model values (space x time)
%   lat_gnss      - [Nx1] vector of mesh latitudes
%   lon_gnss      - [Nx1] vector of mesh longitudes
%   method        - optional argument for interpolation method. Default: 'nearest'
%                 - possible options: 'linear', 'nearest', 'natural'
%
% Outputs:
%   interped_matrix - [MxT] matrix of interpolated model data at GNSS sites over defined GNSS time window

    if nargin < 7
        method = 'nearest';
    end

    num_stations  = length(lon_gnss);
    num_timesteps = size(field_data, 2);
    lon_md = md.mesh.long(:); % Flatten to a 1D vector
    lat_md = md.mesh.lat(:);  % Flatten to a 1D vector
    interped_matrix = NaN(num_stations, num_timesteps);

    for k = 1:num_timesteps
        Z_flat = field_data(:,k);
        interped_matrix(:,k) = griddata(lon_md, lat_md, Z_flat, lon_gnss, lat_gnss, method);
    end
end