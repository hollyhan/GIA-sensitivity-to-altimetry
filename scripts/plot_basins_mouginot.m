function h = plot_basins_mouginot(fpath_basin_boundaries,varargin)

    p = inputParser;
    addParameter(p,'Color',[0.25 0.25 0.25]);
    addParameter(p,'LineWidth',0.8);
    addParameter(p,'Subregions',{'NW','N','NE','CW','CE','SW','SE'});
    parse(p,varargin{:});

    S = shaperead(fpath_basin_boundaries);
    subregions = p.Results.Subregions;
    h = gobjects(numel(subregions),1);

    hold on

    for i = 1:numel(subregions)
        reg = subregions{i};
        ind = strcmp({S.SUBREGION1},reg);
        
        poly_list = polyshape.empty;

        for k = find(ind)
            x = S(k).X;
            y = S(k).Y;

            good = ~isnan(x) & ~isnan(y);
            x = x(good);
            y = y(good);

            if numel(x) < 3
                continue
            end

            poly_list(end+1) = polyshape(x,y,'Simplify',true); %#ok<AGROW>
        end

        if ~isempty(poly_list)
            % 1. Merge all pieces together
            poly_reg = union(poly_list);
            
            % 2. Separate the merged shape into individual distinct pieces (mainland vs islands)
            poly_pieces = regions(poly_reg);
            
            if ~isempty(poly_pieces)
                % Calculate the area of each separate piece
                areas = area(poly_pieces);
                
                % Find the index of the single largest piece (the main ice sheet region)
                [~, largest_idx] = max(areas);
                main_body = poly_pieces(largest_idx);
                
                % 3. Extract coordinates of only this main body
                [xb, yb] = boundary(main_body);
                
                % 4. Plot the clean, island-free boundary
                h(i) = plot(xb, yb, ...
                    'Color', p.Results.Color, ...
                    'LineWidth', p.Results.LineWidth);
            end
        end
    end
end
