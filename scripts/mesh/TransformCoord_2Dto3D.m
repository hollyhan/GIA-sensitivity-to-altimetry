function md_3D = TransformCoord_2Dto3D(md_2D, espg_code)
    % Holly Han (created: Nov. 25th, 2024; Last edited: Nov. 25th, 2024).
    % Convert a 2D spherical coordinate mesh (latitude/longitude) to a 3D surface mesh.
    % This function is essentially same as the function TwoDtoThreeD in the
    % ISSM MATLAB toolbox. But I wrote this version that avoids using 'gdaltransform'
    % and instead directly computes the 3D coordinates from the 2D latitude/longitude.
    %
    % Inputs:   
    %   md_2D: model object with 2D mesh
    %   espg_code: EPSG code of the projection system
    % Outputs:
    %   md_3D: model object with 3D mesh

    % Create a new model object
    md_3D = model();

    % Copy all fields from md_2D to md_3D except the mesh
    fields = fieldnames(md_2D);
    for i = 1:length(fields)
        if ~strcmp(fields{i}, 'mesh')
            md_3D.(fields{i}) = md_2D.(fields{i});
        end
    end

    % Get necessary constants from the provided model
    R = md_2D.solidearth.planetradius;

    % Expected EPSG code for Greenland (Polar Stereographic North)
    expected_epsg = espg_code;

    % Check if the EPSG code in the 2D mesh matches the expected one
    if md_2D.mesh.epsg ~= expected_epsg
        warning(['EPSG code mismatch! Expected EPSG code ' num2str(expected_epsg) ...
                 ' (Polar Stereographic North), but got EPSG: ' num2str(md_2D.mesh.epsg) '.']);

        % assign the right setting for EPSG
        disp(['Setting the EPSG code to ' num2str(expected_epsg) '...']);
        md_2D.mesh.epsg = expected_epsg;
    end

    % Get latitude and longitude from the model
    % check if the mesh is already in the right coordinate system
    if any(~isnan(md_2D.mesh.lat))  && any(~isnan(md_2D.mesh.long)) 
        long = md_2D.mesh.long; % Longitude values of the mesh vertices
        lat = md_2D.mesh.lat;   % Latitude values of the mesh vertices
    else
        % if not, get latitude and longitude from the model
        % Retrieve projection system (typically EPSG 3413 for Greenland)
        proj_id = projcrs(md_2D.mesh.epsg);  % Should be EPSG 3413
    
        % Perform inverse projection to get latitude and longitude
        [lat, long] = projinv(proj_id, md_2D.mesh.x, md_2D.mesh.y);
    end

    % Convert spherical coordinates to Cartesian coordinates
    x = R .* cosd(lat) .* cosd(long); % x-coordinate
    y = R .* cosd(lat) .* sind(long); % y-coordinate
    z = R .* sind(lat);               % z-coordinate

    % Preserve existing mesh connectivity and boundary information
    elements = md_2D.mesh.elements;               % Element connectivity
    edges = md_2D.mesh.edges;
    numberofedges = md_2D.mesh.numberofedges;   
    vertexConnectivity = md_2D.mesh.vertexconnectivity; % Vertex connectivity
    vertexOnBoundary = md_2D.mesh.vertexonboundary;     % Boundary vertex flags

    % Create a new 3D surface mesh
    md_3D.mesh = mesh3dsurface();                  % Initialize a new 3D surface mesh object
    md_3D.mesh.lat = lat;                          % Latitude values
    md_3D.mesh.long = long;                        % Longitude values
    md_3D.mesh.x = x;                              % x-coordinates
    md_3D.mesh.y = y;                              % y-coordinates
    md_3D.mesh.z = z;                              % z-coordinates
    md_3D.mesh.elements = elements;                % Element connectivity
    md_3D.mesh.edges = edges;
    md_3D.mesh.numberofedges = numberofedges;
    md_3D.mesh.numberofelements = length(elements); % Number of elements
    md_3D.mesh.numberofvertices = length(lat);     % Number of vertices
    md_3D.mesh.r = R * ones(md_3D.mesh.numberofvertices, 1); % Radius for all vertices
    md_3D.mesh.area = GetAreas3DTria(md_3D.mesh.elements, ...
                                               md_3D.mesh.x, ...
                                               md_3D.mesh.y, ...
                                               md_3D.mesh.z);
    md_3D.mesh.vertexconnectivity = vertexConnectivity; % Vertex connectivity
    md_3D.mesh.vertexonboundary = vertexOnBoundary;     % Boundary vertex flags
end