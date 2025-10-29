function [md_updated, report] = interpolate_altimetry_to_mesh_massconservative( ...
    h_annual, lat, lon, years, md, X, Y, dhdt_annual, varargin)
% interpolate_altimetry_to_mesh_massconservative
% Resolution-independent, mass-conserving interpolation of altimetry-derived
% ice thickness (h) or rate-of-change (dh/dt) onto an ISSM mesh.
%
% Name–value pairs (optional):
%   'sigma'        - Primary Gaussian kernel width (m) [default = 10e3]
%   'sigma_final'  - Secondary smoother width (m)      [default = 5e3]

tic;
rhoi = 917; % kg/m^3

%% Parse optional args
p = inputParser;
addParameter(p,'sigma',10e3,@(x)isnumeric(x)&&x>0);
addParameter(p,'sigma_final',5e3,@(x)isnumeric(x)&&x>0);
parse(p,varargin{:});
sigma       = p.Results.sigma;
sigma_final = p.Results.sigma_final;

%% Input validation
if nargin < 7
    error('Usage: interpolate_altimetry_to_mesh_massconservative(h_annual,lat,lon,years,md,X,Y,[dhdt_annual])');
end
disp('Running resolution-independent, mass-conserving interpolation...')

if size(lat,1) ~= size(lon,1) || size(lat,2) ~= size(lon,2)
    error('Latitude and longitude dimensions must match.');
end
if size(h_annual,1) == size(lat,2) && size(h_annual,2) == size(lat,1)
    disp('Permuting h_annual dimensions to match lat/lon...');
    h_annual = permute(h_annual,[2 1 3]);
elseif size(h_annual,1) ~= size(lat,1) || size(h_annual,2) ~= size(lat,2)
    error('h_annual spatial dimensions do not match coordinates.');
end

%% Grid info
dX = mean(diff(X));
dY = mean(diff(Y));
cell_area = abs(dX * dY);
fprintf('Projected altimetry grid spacing: %.0f × %.0f m (area %.0f m²)\n', abs(dX), abs(dY), cell_area);
[Xgrid,Ygrid] = meshgrid(X,Y);
cell_area_local = abs(dX * dY) * ones(numel(Xgrid),1);

%% Choose data
if exist('dhdt_annual','var') && ~isempty(dhdt_annual)
    data3D = dhdt_annual; nt = size(dhdt_annual,3); data_label = 'dhdt';
    data3D = double(data3D);
else
    data3D = h_annual;    nt = size(h_annual,3);    data_label = 'h';
    data3D = double(data3D);
end

%% Target mesh
nVerts = md.mesh.numberofvertices;
x = double(md.mesh.x(:));
y = double(md.mesh.y(:));

interp_data     = zeros(nVerts, nt);
mass_src_total  = zeros(nt,1);
mass_dst_total  = zeros(nt,1);

fprintf('Number of mesh vertices: %d \n', nVerts);
fprintf('Smoothing kernels: %.0f km (primary) + %.0f km (final)\n', sigma/1e3, sigma_final/1e3);

%% Per-vertex areas
tri  = md.mesh.elements;
Atri = 0.5 * abs((x(tri(:,2))-x(tri(:,1))) .* (y(tri(:,3))-y(tri(:,1))) - ...
                 (x(tri(:,3))-x(tri(:,1))) .* (y(tri(:,2))-y(tri(:,1))));
Avert = accumarray(tri(:), repmat(Atri/3,3,1), [nVerts 1]);

%% KD-tree setup
kdtree = KDTreeSearcher([x, y]);
[neighbors_primary, distances_primary] = rangesearch(kdtree,[x y],3*sigma);
[neighbors_final,   distances_final]   = rangesearch(kdtree,[x y],3*sigma_final);

%% Main interpolation loop
t_start = tic;
for t = 1:nt
    fprintf('Interpolating timestep %d / %d ...\n',t,nt);
    if t>1
        avg_time = toc(t_start)/(t-1);
        eta_min  = (nt-t+1)*avg_time/60;
        fprintf('  Estimated time remaining: %.2f minutes\n',eta_min);
    end

    slice = data3D(:,:,t);
    mass_src_total(t) = nansum(slice(:))*cell_area*rhoi*1e-12; % Gt
    zv = slice(:); mask = ~isnan(zv);
    F = scatteredInterpolant(Xgrid(mask),Ygrid(mask),zv(mask),'natural','none');
    vals = double(F(x,y));

    % ===== Primary Gaussian smoothing (defensive) =====
    vals_smooth = zeros(nVerts,1);
    for i = 1:nVerts
        neigh = neighbors_primary{i};
        r = distances_primary{i};
        if isempty(neigh)
            vals_smooth(i)=vals(i); continue
        end
        neigh = neigh(:); r = double(r(:));
        neigh = neigh(neigh>=1 & neigh<=nVerts);
        n = min(numel(neigh),numel(r));
        if n==0, vals_smooth(i)=vals(i); continue, end
        neigh = neigh(1:n); r=r(1:n);

        w = exp(-(r.^2)/(2*sigma^2)); w=w(:);
        v = vals(neigh); v=v(:);
        n2 = min(numel(w),numel(v));
        if n2==0, vals_smooth(i)=vals(i); continue, end
        w=w(1:n2); v=v(1:n2);
        s=sum(w); if ~isfinite(s)||s==0, w(:)=1/n2; else, w=w./s; end
        vals_smooth(i)=nansum(w.*v);
    end
    vals = vals_smooth; vals(isnan(vals))=0;
    interp_data(:,t)=vals;
    mass_dst_total(t)=nansum(vals.*Avert)*rhoi*1e-12; % Gt

    % ===== Final Gaussian smoothing (defensive) =====
    vals_final = zeros(nVerts,1);
    for i = 1:nVerts
        neigh = neighbors_final{i};
        r = distances_final{i};
        if isempty(neigh)
            vals_final(i)=vals(i); continue
        end
        neigh = neigh(:); r = double(r(:));
        neigh = neigh(neigh>=1 & neigh<=nVerts);
        n = min(numel(neigh),numel(r));
        if n==0, vals_final(i)=vals(i); continue, end
        neigh = neigh(1:n); r=r(1:n);

        w = exp(-(r.^2)/(2*sigma_final^2)); w=w(:);
        v = vals(neigh); v=v(:);
        n2 = min(numel(w),numel(v));
        if n2==0, vals_final(i)=vals(i); continue, end
        w=w(1:n2); v=v(1:n2);
        s=sum(w); if ~isfinite(s)||s==0, w(:)=1/n2; else, w=w./s; end
        vals_final(i)=nansum(w.*v);
    end

    % Preserve total (area-weighted) mass
    total_before = nansum(vals.*Avert);
    total_after  = nansum(vals_final.*Avert);
    if isfinite(total_after)&&total_after~=0
        vals_final = vals_final*(total_before/total_after);
    end
    interp_data(:,t)=vals_final;
end

%% Local area-weighted renormalization
fprintf('\nApplying local area-weighted renormalization (projected grid)...\n');
Fw=scatteredInterpolant(Xgrid(:),Ygrid(:),cell_area_local,'nearest','nearest');
w_local=Fw(x,y);
for t=1:nt
    vals=interp_data(:,t);
    total_src=nansum(data3D(:,:,t),'all')*cell_area*rhoi*1e-12;
    total_mesh=nansum(vals.*Avert)*rhoi*1e-12;
    if total_mesh>0
        vals=vals.*(total_src./total_mesh).*(w_local./mean(w_local,'omitnan'));
    end
    interp_data(:,t)=vals;
end

%% Mass diagnostics
mass_corr_total=zeros(nt,1);
for t=1:nt
    mass_corr_total(t)=nansum(interp_data(:,t).*Avert)*rhoi*1e-12;
end
report.mass_src_total=mass_src_total;
report.mass_dst_total=mass_dst_total;
report.mass_corr_total=mass_corr_total;
report.mean_error_percent=mean(abs(1-mass_corr_total./mass_src_total)*100);

diff_sum = sum(report.mass_corr_total) - sum(report.mass_src_total);
pct_diff = abs(diff_sum) / abs(sum(report.mass_src_total)) * 100;

fprintf('Global total mass (src vs corr): %.2f vs %.2f Gt — Δ = %.2f Gt (%.2f%%)\n', ...
    sum(report.mass_src_total), sum(report.mass_corr_total), diff_sum, pct_diff);
fprintf('Mean total-mass error after projected-area renormalization: %.2f%%\n',report.mean_error_percent);
fprintf('Mean absolute mismatch: %.2f Gt\n',mean(abs(report.mass_corr_total-report.mass_src_total)));

%% Spatial diagnostic
t_diag=round(nt/2);
slice_src=data3D(:,:,t_diag);
slice_interp=interp_data(:,t_diag);
F_src=scatteredInterpolant(Xgrid(:),Ygrid(:),slice_src(:),'linear','none');
src_on_mesh=F_src(x,y);
mask_valid=~isnan(src_on_mesh)&~isnan(slice_interp);
if sum(mask_valid)>10
    r=corr(src_on_mesh(mask_valid),slice_interp(mask_valid));
    rms_diff=sqrt(mean((src_on_mesh(mask_valid)-slice_interp(mask_valid)).^2));
    fprintf('Diagnostic t=%d: corr=%.3f, RMS=%.2f m\n',t_diag,r,rms_diff);
end
fprintf('Interpolation completed (mass-conservative, resolution-independent).\n');

%% Reconstruct thickness
if strcmp(data_label,'dhdt')
    md.masstransport.spcthickness=zeros(nVerts+1,nt+1);
    fprintf('Initialized spcthickness: %d verts × %d timesteps\n',nVerts,nt+1);
    h_reconstructed=zeros(nVerts,nt+1);
    if ~isempty(h_annual)
        h0=mean(h_annual(:,:,1),'all','omitnan');
    else
        h0=1500;
    end
    fprintf('Using mean initial thickness %.1f m\n',h0);
    h_reconstructed(:,1)=h0;
    for t=1:nt
        h_reconstructed(:,t+1)=h_reconstructed(:,t)+interp_data(:,t);
        h_reconstructed(h_reconstructed(:,t+1)<0,t+1)=0;
    end
    md.masstransport.spcthickness(1:nVerts,:)=h_reconstructed;
else
    md.masstransport.spcthickness=zeros(nVerts+1,nt);
    fprintf('Initialized spcthickness: %d verts × %d timesteps\n',nVerts,nt);
    md.masstransport.spcthickness(1:nVerts,:)=interp_data;
end

md.masstransport.spcthickness(end,:)=years(:)';
disp('Reconstruction of ice thickness on md.masstransport.spcthickness completed.');

md_updated=md;
elapsed_time=toc;
fprintf('Interpolation routine completed in %.2f seconds (%.2f minutes)\n',elapsed_time,elapsed_time/60);
disp('====================================');
end
