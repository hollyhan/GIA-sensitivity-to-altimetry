function [WRMS, mean_err_gnss] = calculate_wrms_and_meanerr(data_gnss, err_gnss, interped_md_data)
% calculate_wrms_and_meanerr computes the weighted RMS error and mean GNSS error.
%
% Inputs:
%   data_gnss          - cell array of GNSS displacement data (not directly used here)
%   err_gnss           - cell array of GNSS errors for each station (same size as data_gnss)
%   interped_md_data   - cell array of model data interpolated to GNSS sites/times
%
% Outputs:
%   WRMS               - weighted root mean square error per station
%   mean_err_gnss      - mean error per station (NaN-omitting)

    num_stations = length(data_gnss);
    WRMS = zeros(num_stations, 1);
    mean_err_gnss = zeros(num_stations, 1);

    for n = 1:num_stations
        sigma = err_gnss{n};  % assumed to match the length of interpolated time series
        w = 1 ./ sigma.^2;

        % WRMS: weighted error between model and GNSS
        disp(['Calculating WRMS for station ', num2str(n)])
        WRMS(n) = sqrt( sum(w .* (data_gnss{n} - interped_md_data{n}).^2) / sum(w) );

        % Mean error, omitting NaNs
        mean_err_gnss(n) = mean(sigma, 'omitnan');
    end
end
