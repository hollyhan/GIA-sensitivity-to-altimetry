function h = plot_basins_rignot(fpath_basin_boundaries,varargin)
    % plot_basins_rignot overlays major Rignot/IMBIE Greenland basins.
    %
    % Usage:
    %   S = shaperead(fpath_basin_boundaries);
    %   plot_basins_rignot(S)
    %
    % Optional:
    %   plot_basins_rignot(S,'Color',[0.35 0.35 0.35],'LineWidth',0.8)

    p = inputParser;
    addParameter(p,'Color',[0.35 0.35 0.35]);
    addParameter(p,'LineWidth',0.8);
    parse(p,varargin{:});

    S = shaperead(fpath_basin_boundaries);

    h = gobjects(length(S),1);

    hold on

    for k = 1:length(S)

        lon = S(k).X;
        lat = S(k).Y;

        % Convert 0-360 lon if needed
        lon(lon > 180) = lon(lon > 180) - 360;

        % Project from lon/lat to polar stereographic north
        [x,y] = ll2psn(lat,lon);

        h(k) = plot(x,y,'-', ...
            'Color',p.Results.Color, ...
            'LineWidth',p.Results.LineWidth);
    end

end