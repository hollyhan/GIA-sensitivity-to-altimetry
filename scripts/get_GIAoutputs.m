function [RSL, thk, B, S, deltaB, deltaS, time_array] = get_GIAoutputs(md)
    % Holly Han (created August 22th, 2024; Last edited in Aug. 22th, 2024).
    % This function takes in a model with transient solutions and
    % postprocess sea-level change relative to the final model timestep.
    % Input: 
    %    - md: model that contains md.results.TransientSolution
    % Output:
    %    - RSL: 2D matrix of relative sea level field through time
    %    - thk: ice thickness 
    %    - deltaB: total change in bed elevation relative to initial time
    %    - deltaS: total change change in sea surface geoid relative to initial time

    % number of timesteps in the ISSM model
    nTime = length(md.results.TransientSolution);
    
    % number of verticies
    nSpace = md.mesh.numberofvertices;

    % Initialize the matrix to store the updated thickness profile
    RSL = zeros(nSpace, nTime);

    % calculate RSL from model
    for j = 1:nTime
        time_array(j) = md.results.TransientSolution(j).time;
        %S(:, j) = md.results.TransientSolution(j).Sealevel; % Elevation of Sea level surface
        S(:, j) = NaN;
        B(:, j) = md.results.TransientSolution(j).Bed; % Elevation of Bedrock
        %thk(:, j) = md.results.TransientSolution(j).Thickness;
        thk(:, j) = NaN;
    end
    deltaB = B - B(:,1);
    %deltaS = S - S(:,1);
    deltaS = NaN;
    %RSL = deltaS - deltaB;
    RSL = NaN;
end