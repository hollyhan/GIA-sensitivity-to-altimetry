function [dice_mass, dice_mass_cum] = get_ice_mass_change(rho, xcoords, ycoords, dhdt, time_array, plot_fig)
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
    % rho - density of material (kg/m^3)
    % plot_fig - (Optional) If true, plot monthly and annual ice mass
    % variations
    %
    % Output:
    % dice_mass - ice mass change per time step (in gigatons)
    % dice_mass_cum - cumulative ice mass change (in gigatons)

    if nargin < 5 || isempty(time_array)
        disp('No time vector provided. Trend will not be calculated')
        time_array = [];
    end

    if nargin < 6
        plot_fig = false;
    end
   
    if length(xcoords) ~= size(dhdt, 1) || length(ycoords) ~= size(dhdt, 2)
        error('The dimensions of dhdt should be [length(xcoords), length(ycoords), length(time)]');
    end
    
    if ~isempty(time_array) && size(dhdt, 3) > 1 && size(dhdt, 3) ~= length(time_array)
        disp(['the length of ''time array'' is: ', num2str(length(time_array))])
        disp(['the length of ''time'' in ''dhdt'' is: ', num2str(size(dhdt, 3))])
        error('The length of the input ''time_array'' and length of ''time'' in ''dhdt'' must be the same.')
    end

    dx = abs(diff(xcoords(1:2)));
    dy = abs(diff(ycoords(1:2)));
    cell_area = double(dx * dy);

    disp('ice altimetry original data resolution:')
    disp(['dx (m): ', num2str(dx)]);
    disp(['dy (m): ', num2str(dy)]); 
    disp(['cell_area (m^2): ', num2str(cell_area)]);

    volume_change = dhdt * cell_area;
    mass_change = volume_change * rho * 1.0e-12;

    dice_mass = squeeze(sum(mass_change, [1, 2], 'omitnan'));

    % Cumulative mass change time series
    if size(dhdt, 3) > 1
        dice_mass_cum = cumsum(dice_mass, 'omitnan');
    else
        dice_mass_cum = dice_mass;
    end

    total_mass_change = dice_mass_cum(end);
    mean_rate = mean(dice_mass, 'omitnan');

    if ~isempty(time_array)
        fprintf('Total cumulative mass change: %.1f Gt (%d–%d)\n', ...
            total_mass_change, round(time_array(1)), round(time_array(end)));
    else
        fprintf('Total cumulative mass change: %.1f Gt\n', total_mass_change);
    end

    fprintf('Mean mass change per time step: %.1f Gt\n', mean_rate);

    if length(time_array) > 1
        p = polyfit(time_array, dice_mass, 1);
        trend_per_time = p(1);
        intercept = p(2);
        disp(['Trend in mass change per time step: ', num2str(trend_per_time), ' Gt/unit time']);
        dmass_fit = polyval(p, time_array);

        if plot_fig
            figure;
            plot(time_array, dice_mass, 'o', 'DisplayName', 'Mass Change Data');
            hold on;
            plot(time_array, dmass_fit, '-r', 'DisplayName', ['Linear Fit: ', num2str(trend_per_time, '%.2f'), ' Gt/unit time']);
            xlabel('Time');
            ylabel('Mass Change (Gt)');
            title('Linear Trend in Mass Change Per Time Step');
            legend('show');
            grid on;
        end
    end
end