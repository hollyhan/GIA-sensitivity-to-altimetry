% settings_gia_parameterization.m
% Holly Han (created Oct. 14th, 2024; Last edited on Nov. 14th, 2024).
% This script contains settings for parameterizing the solid earth module
% The options defined here will be used to parameterize model in the code
% 'run_gia.m' Any extra parameterization done in this script other than
% already defined will be used to override the default values of
% corresponding pameterizations done in 'run_gia.m'

% Whether to use the whole domain or GNSS sites only for GIA calculation
% If true, the GIA calculation will be done for the whole domain using the ISSM
% If false, the GIA calculation will be done for the GNSS sites only using the
% GIA Greens Function
use_whole_domain = false;

% Whether Love Numbers are provided from an external file
% If false, LNs need to be calculated in Step 2 in 'main_script.m'
% If true, the file will be read from the path defined below.
earthmodel_external = false;
fpath_earthmodel = '../data_model/earth_structure/love_numbers/lovenumbers_VE_incompressible_Maxwell_LT120_umv0.03_lmv2.mat';
rheology_choice = 'Maxwell'; % Choice of rheological model in String (e.g. 'Maxwell', 'EBM', 'SimpleBM','Elastic'). This variable will be used in the model file name.

% Solution parameters
% Enable Gravity, Rotation and Deformation of the Solid-Earth to be
% computed
mindeg = 1;
degmax_sh = 1e4; % spherical harmonics max degree to solve for

% Sea-level solver settings
compute_GRD_patterns = true; % for 'md.solidearth.settings.isgrd'
enable_selfattraction = true; % for 'md.solidearth.settings.selfattraction'
enable_rotation = false; % for 'md.solidearth.settings.rotation'

% Deformational effects settings. Separately done for each period of interest.
enable_elastic_deformation_present = true; % for 'md.solidearth.settings.elastic' for the present-day run setup
enable_viscous_deformation_present = true; % for 'md.solidearth.settings.viscous'

% Set up simulation interval
time_interval_present = 1; % for 'md.timestepping.time_step'

% Requested outputs for model solve. If empty, only default outputs
% will be written out. 
requested_outputs_solidearth = {'Sealevel','Bed','SealevelBarystaticIceLoad','SealevelBarystaticIceArea','SealevelGRD'};
requested_outputs_masstransport = {'Thickness','DeltaIceThickness'};

% Whether to write out the solved model object into a file
saveModel = true;% Logical variable to decide whether to write the love numbers to a file
fpath_output = '../results/model_objects_saved/'; % Path to which love numbers will be saved if saveModel is true
