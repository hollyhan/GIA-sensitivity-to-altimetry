%% MATLAB script
function [md_global, md_world, md_regional_2D, transfer] = createGlobalMesh(md_regional)
    % createGlobalMesh.m
    % Creates a global mesh for GIA calculations by creating a world mesh
    % (excluding the regional domain) and merging it with the regional mesh.
    %
    % Inputs:
    %   md_regional - Regional mesh model for Greenland
    %
    % Outputs:
    %    - md_global - Final 3D global mesh model ready for GIA calculations
    %    - md_world - 3D global mesh where regional mesh provided by md_regional is extruded
    %    - md_regional_2D: model object with 2D regional mesh
    %    - transfer: the list of indices to refer to the index of regional mesh on the merged global mesh, size of md_global_extruded_3D.mesh.numberofverticies+md_regional_3D.mesh.numberofvertices
	% col 1 is for unmerged mesh, col2 is merged mesh. output of the function 'removeduplicatevertices'
    
    % First, check if the regional mesh is 2D or 3D. If in 3D, project onto a 2D surface
    switch class(md_regional.mesh)
        case 'mesh3dsurface' % This is not supported yet
            % Extract the longitude and latitude from 3D coordinates
            % This step needs to be implemented based on your coordinate system
            md_regional_3D = md_regional;
        case 'mesh2d'
            md_regional_2D = md_regional;
            disp('Transforming 2D regional mesh to 3D')
            md_regional_3D = TransformCoord_2Dto3D(md_regional_2D, 3413);
    end

    % Step 1: Create world mesh (excluding regional domain)
    disp('Creating world mesh (excluding regional domain)...');
    md_world = createWorldMesh(md_regional_2D);

    % Step 2: Merge regional and world meshes
    disp('Merging regional and world meshes...');
    [md_global, transfer] = mergeMeshes(md_regional_3D, md_world);

    % Transfer necessary fields from regional to global mesh
    disp('Transferring fields from regional to global mesh...');
    md_global = transfer_fields_to_global(md_regional_2D, md_world, md_global, 'geometry.thickness', transfer, 'regional');
    md_global = transfer_fields_to_global(md_regional_2D, md_world, md_global, 'geometry.surface', transfer, 'regional');
    md_global = transfer_fields_to_global(md_regional_2D, md_world, md_global, 'geometry.base', transfer, 'regional');
    md_global = transfer_fields_to_global(md_regional_2D, md_world, md_global, 'geometry.bed', transfer, 'regional');
    md_global = transfer_fields_to_global(md_regional_2D, md_world, md_global, 'mask.ice_levelset', transfer, 'regional');
    md_global = transfer_fields_to_global(md_regional_2D, md_world, md_global, 'mask.ocean_levelset', transfer, 'regional');
    md_global = transfer_fields_to_global(md_regional_2D, md_world, md_global, 'masstransport.spcthickness', transfer, 'regional');

    disp('Global mesh creation completed successfully');
end 