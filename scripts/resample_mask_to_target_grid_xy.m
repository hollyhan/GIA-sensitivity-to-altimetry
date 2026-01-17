function mask_resampled = resample_mask_to_target_grid_xy(mask_data, x_mask, y_mask, x_tgt, y_tgt)
    % x_mask, y_mask: [Ny_mask x Nx_mask] in meters (EPSG:3413)
    % x_tgt,  y_tgt : [Ny_tgt  x Nx_tgt ] in meters (EPSG:3413) for the dh/dt grid
    % mask_data     : [Ny_mask x Nx_mask] or [Ny_mask x Nx_mask x Nt], logical/0-1

    % infer 1D target vectors (assumes rectilinear grid in EPSG:3413)
    % take first row for x, first column for y
    xvec = x_tgt(1, :);   % 1 x Nx_tgt, monotonic
    yvec = y_tgt(:, 1);   % Ny_tgt x 1,  monotonic

    % enforce monotonic increasing (discretize requires this)
    flipX = false; flipY = false;
    if numel(xvec) > 1 && xvec(2) < xvec(1)
        xvec  = fliplr(xvec);    % now increasing
        flipX = true;
    end
    if numel(yvec) > 1 && yvec(2) < yvec(1)
        yvec  = flipud(yvec);    % now increasing
        flipY = true;
    end

    % build pixel edges (half-cell padding at ends)
    dx = median(diff(xvec));  dy = median(diff(yvec));
    x_edges = [xvec - dx/2, xvec(end) + dx/2]; % NOTE the comma (horzcat)
    y_edges = [yvec - dy/2; yvec(end) + dy/2]; % NOTE the semicolon (vertcat)

    % flatten mask coords
    xm = x_mask(:);
    ym = y_mask(:);

    % bin indices for each mask pixel center
    ix = discretize(xm, x_edges);
    iy = discretize(ym, y_edges);

    % if we flipped the axes earlier, remap the indices back
    if flipX
        ix = numel(xvec) - ix + 1;
    end
    if flipY
        iy = numel(yvec) - iy + 1;
    end

    in = ~isnan(ix) & ~isnan(iy);
    lin_idx = sub2ind([numel(yvec), numel(xvec)], iy(in), ix(in));

    Ny = numel(yvec); Nx = numel(xvec);
    if ndims(mask_data) == 3
        Nt = size(mask_data, 3);
        mask_resampled = false(Ny, Nx, Nt);
        for t = 1:Nt
            v = mask_data(:,:,t);
            v = v(:);
            v = v(in) > 0;  % logical
            % counts per target cell
            ice_counts   = accumarray(lin_idx, double(v), [Ny*Nx, 1], @sum, 0);
            total_counts = accumarray(lin_idx, 1,         [Ny*Nx, 1], @sum, 0);
            frac = ice_counts ./ max(total_counts,1);
            mask_resampled(:,:,t) = reshape(frac > 0.5, Ny, Nx);
        end
    else
        v = mask_data(:);
        v = v(in) > 0;
        ice_counts   = accumarray(lin_idx, double(v), [Ny*Nx, 1], @sum, 0);
        total_counts = accumarray(lin_idx, 1,         [Ny*Nx, 1], @sum, 0);
        frac = ice_counts ./ max(total_counts,1);
        mask_resampled = reshape(frac > 0.5, Ny, Nx);
    end
end
