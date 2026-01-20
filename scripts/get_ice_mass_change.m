function [dice_mass, dice_mass_cum] = get_ice_mass_change(rhoi, xcoords, ycoords, dhdt, time_array, plot_fig)
    % get_ice_mass_change.m
    % Holly Han (created: Nov. 5th, 2024; Last edited: Nov. 5th, 2024).
    % Takes monthly ice elevation change field on a polarstereographic
    % projection and calculates total ice mass changes over the domain.
    % This function assumes a constant grid spacing in x and y direction.
    %
    % Input:
    % xcoords - X coordinates of a domain (e.g. polarstereograhic projection)
    % ycoords - Y coordinates of a domain
    % dhdt - ice elevation changes on the domain (in meters) in size of
    %        [length(xcoords) length(ycoords) length(time)]
    % time_array (Optional) - time vector that corresponds to timestamps of 'dhdt'
    % rhoi - density of ice (km/m^3)
    % plot_fig - (Optional) If true, plot monthly and annual ice mass
    % variations
    %
    % Output:
    % dice_mass - ice mass change (in gigatons)
    % dice_mass_cum - cumulative ice mass change (in gigatons)

    % Check if time argument is provided
    if nargin < 5 || isempty(time_array)
        % If not provided, use a default time vector based on data size (assuming monthly data)
        disp('No time vector provided. Trend will not be calculated')
        time_array = [];
    end

    if nargin < 6
        plot_fig = false;  % Default plot flag to false if not provided
    end
   
    % Check the data size
    if length(xcoords) ~= size(dhdt, 1) || length(ycoords) ~= size(dhdt, 2)
        error('The dimensions of dhdt should be [length(xcoords), length(ycoords), length(time)]');
    end
    
    if ~isempty(time_array) && size(dhdt, 3) > 1 && size(dhdt, 3) ~= length(time_array)
        disp(['the length of ''time array'' is: ', num2str(length(time_array))])
        disp(['the length of ''time'' in ''dhdt'' is: ', num2str(size(dhdt, 3))])
        error('The length of the input ''time_array'' and length of ''time'' in ''dhdt'' must be the same.')
    end

    % Calculate grid cell area (assuming polar stereographic grid)
    dx = abs(diff(xcoords(1:2)));
    dy = abs(diff(ycoords(1:2)));  % grid spacing in y-direction
    cell_area = double(dx*dy);  % Area of each cell in m^2

    disp('ice altimetry original data resolution:')
    disp(['dx (m): ', num2str(dx)]);
    disp(['dy (m): ', num2str(dy)]); 
    disp(['cell_area (m^2): ', num2str(cell_area)]);

    % Calculate volume change per cell
    volume_change = dhdt * cell_area;  % Volume change per cell, with same dimensions as dh (m^3)

    % Convert volume change to mass change
    mass_change = volume_change * rhoi * 1.0e-12 ;  % Mass change per cell (Gt)

    % Sum mass change over all cells for each time step
    dice_mass = squeeze(sum(mass_change, [1, 2], 'omitnan')); % Sum over latitude and longitude

    % Cumulative mass change time series
    if length(time_array) > 1
        dt = diff(time_array);
        % Handle uneven time spacing (e.g., 1992, 1993.5, 1995)
        mean_dt = [dt(1); dt(:)];  % simple alignment
        dice_mass_cum = cumsum(dice_mass .* mean_dt, 'omitnan');
        dice_mass_cum = dice_mass_cum - dice_mass_cum(1); % relative to start
    else
        dice_mass_cum = 0;
    end

    total_mass_change = dice_mass_cum(end);
    mean_rate = mean(dice_mass, 'omitnan');
    fprintf('Total cumulative mass change: %.1f Gt (%d–%d)\n', total_mass_change, round(time_array(1)), round(time_array(end)));
    fprintf('Mean mass change rate: %.1f Gt/yr\n', mean_rate);

    % Perform linear regression to get trend if there is more than one time element
    if length(time_array) > 1
        % Fit a linear trend
        p = polyfit(time_array, dice_mass, 1);
        % Extract the slope (Gt/yr) and intercept
        trend_per_time = p(1);  % Slope, representing the trend in Gt/yr
        intercept = p(2);       % Y-intercept of the fitted line
        disp(['Trend (mass change rate): ', num2str(trend_per_time), ' Gt/yr']);
        dmass_fit = polyval(p, time_array);  % Calculate fitted values using the trend line

        if plot_fig
            figure;
            plot(time_array, dice_mass, 'o', 'DisplayName', 'Mass Change Data');  % Original data points
            hold on;
            plot(time_array, dmass_fit, '-r', 'DisplayName', ['Linear Fit: ', num2str(trend_per_time, '%.2f'), ' Gt/unit time']);  % Fitted line
            xlabel('Time');
            ylabel('Mass Change (Gt)');
            title('Linear Trend in Rate of Mass Change');
            legend('show');
            grid on;
        end
    end
end