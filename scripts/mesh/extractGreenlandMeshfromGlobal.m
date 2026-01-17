function [meshG, indexG] = extractGreenlandMeshfromGlobal(md)
% extractGreenlandMeshfromGlobal
% Extracts Greenland-only portion of a global ISSM mesh.
% Note: this script requires Arctic Mapping tools
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. (2017). Antarctic Mapping Tools for Matlab. 
%         Computers & Geosciences, 104, 151–157. Elsevier BV. https://doi.org/10.1016/j.cageo.2016.08.003
% Returns a struct with fields:
%   elements   (Nx3)
%   lat, long   (vertex coordinates)
%   x, y (projected coords in AMT polar stereographic)

    % --- geographic coordinates ---
    lat = md.mesh.lat;
    lon = md.mesh.long;

    % ensure lon is in [-180,180]
    lon(lon > 180) = lon(lon > 180) - 360;

    % --- Greenland bounding box ---
    inG = (lat > 58 & lat < 85 & lon > -80 & lon < -10);

    % --- elements whose *all 3* nodes lie in Greenland ---
    elem = md.mesh.elements;
    elem_inG = all(inG(elem), 2);
    elementsG_old = elem(elem_inG, :);

    % --- reindex vertices old → new ---
    old2new = zeros(size(inG));
    old2new(inG) = 1:nnz(inG);
    elementsG = old2new(elementsG_old);

    % --- subset fields ---
    latG  = lat(inG);
    lonG  = lon(inG);

    % --- project to AMT polar stereographic ---
    [xpsnG, ypsnG] = ll2psn(latG, lonG);

    % --- assemble output struct ---
    meshG.elements = elementsG;
    meshG.lat = latG;
    meshG.long = lonG;
    meshG.x = xpsnG;
    meshG.y = ypsnG;
    meshG.numberofvertices = nnz(inG);
    meshG.numberofelements = size(elementsG,1);
    indexG = inG;

end
