function md_solved = run_gia(md_global, ht, kt, lt, tht, tkt, tlt, pmtf1, pmtf2, time, output_fname)
    % Holly Han (created Oct. 14th, 2024; Last edited on Aug. 6th, 2025).
    % This function parameterizes the global model and runs the GIA solver
    % given a prescribed ice loadng history.
    % Inputs:
    %    - md_global: model object that contains a global mesh in md.mesh
    %    - ht: load Love number for radial displacement
    %    - kt: load Love number for gravitational potential perturbation
    %    - lt: load Love number for horizontal displacements
    %    - tht: tidal load Love number (deg 2)
    %    - tkt: tidal load Love number (deg 2)
    %    - tlt: tidal load Love number (deg 2)
    %    - pmtf1: Love number for degree 2 (secular polar motion)
    %    - pmtf2: Love number for chandler wobble
    %    - time: time vector at which Love numbers are solved for
    %    - output_fname: name of the output file [optional]
    % Output:
    %    - md_solved: model object solved for GIA given loading hiscoty

    disp("Calling the script 'run_gia'.m")
    md = md_global;
    md.miscellaneous.name = sprintf('md_global_gia.mat');

    % The below fields need to technically exist and be the right size for
    % to enable masstransport
    md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
    md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
    md.initialization.vx=zeros(md.mesh.numberofvertices,1);
    md.initialization.vy=zeros(md.mesh.numberofvertices,1);
    md.initialization.sealevel = zeros(md.mesh.numberofvertices,1);
    md.smb.mass_balance=zeros(md.mesh.numberofvertices,1);
    md.masstransport.min_thickness = 1e-16;

    % Computation flags
    md.cluster=generic('name',oshostname(),'np',3);

    % Initialize some fields before reading in values from config options
    md.solidearth.settings.elastic = 0; % enable elastic deformation from surface loading
    md.solidearth.settings.viscous = 0; % enable viscous deformation from surface loading
    md.solidearth.settings.isgrd = 0; % compute GRD patterns
    md.solidearth.settings.selfattraction = 0; % enable surface mass load to perturb the gravity field
    md.solidearth.settings.rotation = 0; % enable polar motion to feedback on the GRD field


    % Rest of the parameterization of the sea-level calculation by default
    md.solidearth.settings.grdmodel=1;
    md.solidearth.settings.reltol=1e-3;
    md.solidearth.settings.abstol=NaN;
    md.solidearth.settings.sealevelloading=1;
    md.solidearth.settings.ocean_area_scaling=0;

    % To run the masstransport and sea-level change cores only
    md.transient.issmb=0;
    md.transient.isstressbalance=0;
    md.transient.isthermal=0;
    md.transient.ismasstransport=1;
    md.transient.isslc=1;
    
    % The settings above will be overriden if defined in
    % 'settings_gia_parameterization.m'
    % Read in parameter settings from config options
    disp('Reading in the parameterization settings for the run');
    run('settings_gia_parameterization.m');

    % Assign values from the config options
    md.timestepping.start_time = md.masstransport.spcthickness(end,1);
    md.timestepping.time_step = md.masstransport.spcthickness(end,2) - md.masstransport.spcthickness(end,1);
    md.timestepping.final_time = md.masstransport.spcthickness(end,end);
    if enable_elastic_deformation_present
        disp('Elastic deformation is enabled')
        md.solidearth.settings.elastic = 1;
    end
    if enable_viscous_deformation_present
        disp('Viscous deformation is enabled')
        md.solidearth.settings.viscous = 1;
    end

    % Assign config options for the sea-level solver settings
    if compute_GRD_patterns
        disp('GRD computation is enabled')
        md.solidearth.settings.isgrd = 1;
    end

    if enable_selfattraction
        disp('Self-attraction is enabled')
        md.solidearth.settings.selfattraction = 1;
    end

    if enable_rotation
        disp('Rotation is enabled');
        md.solidearth.settings.rotation = 1;
    end

    % Love Numbers parameterization
    mindeg = 1;
    maxdeg = degmax_sh;
    % lt = [];
    if (md.solidearth.settings.viscous == 0 && md.solidearth.settings.elastic == 1)
        % Alternative option to get the elastic loading from love numbers:
        % md.solidearth.lovenumbers=lovenumbers('maxdeg',degmax_sh)
        % But here we are using the numbers calculated from
        % 'calculate_lovenumbers.m'
        md.solidearth.lovenumbers.h = ht(:,1);
        md.solidearth.lovenumbers.k = kt(:,1);
        md.solidearth.lovenumbers.l = lt(:,1);
        md.solidearth.lovenumbers.th = tht(1:maxdeg+1,1);
        md.solidearth.lovenumbers.tk = tkt(1:maxdeg+1,1);
        md.solidearth.lovenumbers.tl = tlt(1:maxdeg+1,1);
        md.solidearth.lovenumbers.pmtf_colinear = pmtf1(:,1);
        md.solidearth.lovenumbers.pmtf_ortho = pmtf2(:,1);
        md.solidearth.lovenumbers.timefreq = 0;
    elseif (md.solidearth.settings.viscous == 1 && md.solidearth.settings.elastic == 1)
        md.solidearth.lovenumbers.h = ht;
        md.solidearth.lovenumbers.h(1:mindeg,:) = 0;
        md.solidearth.lovenumbers.k = kt;
        md.solidearth.lovenumbers.k(1:mindeg,:) = -1;
        md.solidearth.lovenumbers.l = lt;
        md.solidearth.lovenumbers.l(1:mindeg,:) = 0;
        md.solidearth.lovenumbers.th = tht(1:maxdeg+1,:);
        md.solidearth.lovenumbers.tk = tkt(1:maxdeg+1,:);
        md.solidearth.lovenumbers.tl = tlt(1:maxdeg+1,:);
        md.solidearth.lovenumbers.pmtf_colinear = pmtf1;
        md.solidearth.lovenumbers.pmtf_ortho = pmtf2;
        md.solidearth.lovenumbers.timefreq = time;
    end

    md.solidearth.requested_outputs = requested_outputs_solidearth;
    md.masstransport.requested_outputs = requested_outputs_masstransport;

    % solve model
    disp('Solving for GIA ...');
    md_solved = solve(md,'Transient');

    if saveModel
        if isempty(output_fname)
            output_fname = 'md_gia_solved.mat';
        else
            fname = sprintf('%s.mat', output_fname);
        end
     
        fpath_full = fullfile(fpath_output, fname);
        
        % Create directory if it doesn't exist
        if ~exist(fpath_output, 'dir')
            mkdir(fpath_output);
            disp(['Created output directory: ', fpath_output]);
        end
        
        save(fpath_full, 'md_solved');
        disp(['Solved model object saved to a file: ', fpath_full]);
    else
        disp('saveModel is disabled. Not writing into a file.')
    end

end