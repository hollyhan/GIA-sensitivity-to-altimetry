function md_rest_of_world = createWorldMesh(md_greenland_2D)
    % createWorldMesh.m
    % Holly Han (created Nov. 27th, 2024; Last edited on Nov. 27th, 2024).
    %
    % This function generates a mesh for the rest of the world.
    % 
    % Inputs:
    %    - md_greenland_2D: model object with 2D Greenland mesh
    % Output:
    %    - md_rest_of_world: model object with 3D mesh for the rest of the world

    % Display message
    disp('Generating mesh for the rest of the world...');
    
    % Extract boundary points from the Greenland mesh
    bdry_points = [md_greenland_2D.mesh.segments(:,1); md_greenland_2D.mesh.segments(1,1)];
    
    % Get latitude and longitude of the boundary points
    dom.lat = md_greenland_2D.mesh.lat(bdry_points);
    dom.long = md_greenland_2D.mesh.long(bdry_points);

    % Create the rest of the world mesh using meshing_global_extruded
    coastres=300e3; % resolution at present-day coastline (in meters)
    land_res=2000e3; % resolution inland (in meters)
    ocean_res=4000e3; % resolution in the ocean (in meters)
    smoothness=0; % how smoothly we transition between resolutions; '2' might be a sweet spot
    polecapS=32; % defines how far south we use polar projection to mesh Antarctica
    boundrefine=2; % how many iterations we do in the meshing to make sure different partitions transition smoothly into one another
    smallconts=100; % reduce coastline contours with initially this many points or less to just triangles (mostly done to avoid lots of points around islands)
    minpoint=300; % discard contours with this many points or less (only deal with bigger land masses)
    %Refine.H=md_greenland_2D.mesh.sHi/117;
    %Refine.H(md_greenland_2D.mesh.sHi==0)=NaN;
    %Refine.lat=md_greenland_2D.mesh.lat;
    %Refine.long=md_greenland_2D.mesh.long;
    %Refine.H=(md_greenland_2D.mesh.lat-60).^1.5.*(md_greenland_2D.mesh.lat>60);Refine.H(Refine.H==0)=NaN;
    %Refineerr=5;

    %'arctic_divide' format:[long;lat], adjust if the mesher crashes
    %telling you that the basin domain containing the Americas might
    %intersect with coastlines; identify points of the boundary by their
    %latitude
    Nlim=[-180 -124.88 -114.08 -112.21  -108.8 -106.15 -102.61 -99.082 -96.082-.5 -95.682 -70.95 -60.622 -49.891 -2.2996  17.563  56.223 74.153 108.88  135.07 168.67 180
            72  70.674  68.874  68.408   68.7   68.98   68.68   69.474  72.073+.5  74.206  74.206 67.275  52.814  64.211  74.14   76.541 80.671 82.805  78.14  72.473 72];

    md_rest_of_world_3D = world_mesher('coastres', coastres, 'land_res', land_res, ...
        'ocean_res', ocean_res, 'smoothness', smoothness, ...
        'polecapS', polecapS, 'arctic_divide', Nlim, ...
        'boundrefine', boundrefine, 'smallconts', smallconts, ...
        'minpoint', minpoint, 'extradom', dom, 'plotting', 0);
        
    % Assign the generated mesh to the output variable
    md_rest_of_world = md_rest_of_world_3D;
    
    % Plot for debugging
    %plotmodel(md_rest_of_world,'data', md_rest_of_world.mask.ocean_levelset,'edgecolor','k','coord','xy');caxis([-1 1])
    %set(gca,'clipping','off') 
end
