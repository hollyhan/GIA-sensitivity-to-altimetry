function md_destination = transfer_fields_to_global_mesh(md_source_regional, md_source_global_extruded, md_destination, fieldname, ind, option)
% transfer_fields_to_global_mesh.m
% Holly Han (created Nov. 27th, 2024; Last edited on April 2nd, 2025).
%
% This function transfers field data from a 3D globally extruded mesh model
% and a regional mesh model to a 3D merged global mesh model. 
%
% Parameters:
%   md_source_regional (struct): The source regional mesh model containing
%       the original field data to be transferred. Mesh can be 2D (recommended) or 3D.
%       Important note: 3D meshes created by mesh3dsurface() will not work for the function project2d
%       because there is no information on 'numberoflayers' and 'numberofverticies2d'
%   md_source_global_extruded (struct): The source mesh global-extruded model containing
%        the original field data to be transferred. Needs to be in 3D.
%   md_destination (struct): The destination mesh model (i.e., global-merged mesh)
%       onto which data fields are transferred from the source models.
%   field_name (string): The source field data to be transferred. The format needs to include
%       'field' of the model instance and 'subfield' of the field where md.(field).(subfield) 
%       e.g., 'geometry.thickess'.
%   ind (array): An index mapping array that specifies how vertices
%                from the regional mesh correspond to vertices on the global mesh.
%                i.e. 'transfer' array outputted from function 'mergeMeshes'
%   option (string): The mesh from which the field data is being transferred.
%       Available options are 'regional', 'global_extruded', or 'both' ('both is not available for now.').
%
% Returns:
%   md_global_merged (struct): The updated global mesh model with transferred field data/
%
% Example:
%   md_global_merged = transfer_fields_to_global(md_regional, md_global_extruded, md_global_merged, geometry.thickness, ind, 'regional');
%
% Usage:
%   This function is typically used in mesh processing workflows where
%   field data needs to be consolidated from multiple mesh sources into
%   a single global mesh model.

% Split the field name into parts
parts = strsplit(fieldname, '.');

% Field to be interpolated
var_tmp = md_source_regional.(parts{1}).(parts{2});

% Get number of vertices of the source regional mesh in 2D
switch class(md_source_regional.mesh)
case 'mesh3dsurface'
    if ~isprop(md_source_regional.mesh, 'numberofvertices2d')
        error("'mesh.numberofvertices2d' does not exist. Cannot project 3D onto 2D. Use the 3D regional mesh.");
    else
        % Project the regional source model mesh onto 2D if 3D.
        % Special treatment for spcthickness because the last index is a time array
        if strcmp(parts{1}, 'masstransport') && strcmp(parts{2}, 'spcthickness')
            if length(var_tmp) > n_vertices_regional2D
                var_source_regional = project2d(md_source_regional, var_tmp(1:end-1,:), 1);
            else
                var_source_regional = var_tmp;
            end
        else
            var_source_regional = project2d(md_source_regional, var_tmp, 1);
        end
    end
case 'mesh2d'
    n_vertices_regional2D = md_source_regional.mesh.numberofvertices;
    % Special treatment for spcthickness because the last index is a time array
    if strcmp(parts{1}, 'masstransport') && strcmp(parts{2}, 'spcthickness')
        if length(var_tmp) > n_vertices_regional2D
            var_source_regional = var_tmp(1:end-1,:);
        else
            var_source_regional = var_tmp;
        end
    else
        var_source_regional = var_tmp;
    end
end

% Get the number of vertices and timestamps for the variable
[n_vertices_regional2D, nt_regional] = size(var_source_regional);

if strcmp(option, 'regional')
    if strcmp(fieldname, 'masstransport.spcthickness')
        md_destination.(parts{1}).(parts{2}) = zeros(md_destination.mesh.numberofvertices + 1, length(md_source_regional.(parts{1}).(parts{2})(1,:))); % Use NaN for initialization
        md_destination.(parts{1}).(parts{2})(end, :) = md_source_regional.(parts{1}).(parts{2})(end, :);
    end
    % Transfer fields from the regional mesh to the global mesh
    for t = 1:nt_regional
        for i = 1:n_vertices_regional2D
            md_destination.(parts{1}).(parts{2})(ind(i + md_source_global_extruded.mesh.numberofvertices, 2), t) = md_source_regional.(parts{1}).(parts{2})(i, t);
        end
    end
    disp(['Field "', fieldname, '" transferred from the regional mesh successfully!']);
end

% Transfer values from the global extruded mesh to the destination mesh
if strcmp(option, 'global_extruded')
    [n_vertices, nt_global] = size(md_source_global_extruded.(parts{1}).(parts{2}));
    % Transfer fields from the extruded global mesh to the global mesh
    for t = 1:nt_global
        for i = 1:md_source_global_extruded.mesh.numberofvertices
            md_destination.(parts{1}).(parts{2})(ind(i, 2), t) = md_source_global_extruded.(parts{1}).(parts{2})(ind(i, 1), t);
        end
    end
    disp(['Field "', fieldname, '" transferred from the extruded global mesh successfully!']);
end

if strcmp(option, 'both') % This is currently unavailable.
    [n_vertices, nt_global] = size(md_source_global_extruded.(parts{1}).(parts{2}));
    % Transfer fields from both the extruded global mesh and the regional mesh to the global mesh
    for t = 1:nt_global
        for i = 1:md_source_global_extruded.mesh.numberofvertices
            md_destination.(parts{1}).(parts{2})(ind(i, 2), t) = md_source_global_extruded.(parts{1}).(parts{2})(ind(i, 1), t);
        end
    end
    for t = 1:nt_regional
        for i = 1:n_vertices_regional2D
            md_destination.(parts{1}).(parts{2})(ind(i + md_source_global_extruded.mesh.numberofvertices, 2), t) = md_source_regional.(parts{1}).(parts{2})(i, t);
        end
    end
    disp(['Field "', fieldname, '" transferred from both the extruded global mesh and the regional mesh successfully!']);
end
