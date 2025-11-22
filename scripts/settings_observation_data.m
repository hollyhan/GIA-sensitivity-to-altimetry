% settings_observation_data.m
% Holly Han (created Oct. 15th, 2024; Last edited on Nov. 14th, 2024).
% Note: the processed data will cover different time period depending on
% the source: 'measureItsLive' will cover 1993-2023, and 'DTU' will cover
% 2003-2022.

% Set up some constants
rhoi = 917.0;  % ice density in kg/m3. Used in Nilson et al. Also check consistency across other models
rhoo = 1000.0; % ocean density in kg/m3
r_earth = 6371000.0; % radius of the earth typically used in ISSM

% path to love numbers
fpath_love_numbers = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/lovenumbers/lovenumbers_VE_compressible_Maxwell_LT70_umv1_lmv1.mat';

% path to mesh model
fpath_mesh_model_regional = '/Users/kyhan/Desktop/LIA-Project/Pre-process/issm_model/source_file_josh/Param_UW.mat';
fpath_mesh_model_global = '/Users/kyhan/Desktop/LIA-Project/Pre-process/IceHistory/md_global_ready_for_solve.mat';
fpath_mesh_model_regional_refined = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results/mesh';
fpath_results_general = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results';
fpath_results_figures = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/results/figures';

%--------------------- Altimetry Data Setting --------------------%
% Name of the ice elevation data to use. Possible options are, (1)
% (1)'DTU', (2)'measureItsLive'
data_name_altimetry = 'measureItsLive';

% If data_name is set to 'measureItsLive', also need to choose a firn model
% Available options are, (1) 'GEMB' and (2) 'GSFC'
firnmodel_name = 'GEMB';

% Settings for the geodetic observations for ice thickness change for the
% contemporary period
% Notes on data:
% Khan et al provides "mean monthly elevation change rates corrected for GIA; elastick uplift and firn compaction in water equivalent"
fpath_ice_elev_DTU = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/khan_et_al/Greenland_dhdt_mass_1kmgrid.nc' ;% Khan et al
fpath_iceErr_DTU = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/khan_et_al/Greenland_dhdt_mass_err_1kmgrid.nc';

fpath_ice_elev_measureItsLive = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/nilsson_et_al/Synthesis_GrIS_Perf_r1920m_1992_2023_ESSD_Compliant_v2.nc';
fpath_fac_gemb = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/nilsson_et_al/GEMB_FAC_SMB_ALT_GRID.h5'; % FAC model Glacier Energy and Mass Balance FDM (GEMB) version 1.2 (Gardner et al., 2023) 
fpath_fac_gsfc = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/altimetry/nilsson_et_al/GSFC_GrIS_FAC_SMB_ALT_GRID.h5';% FAC model Goddard Space Flight Center FDM version 1.2.1 (Medley et al., 2022).


%----------------------- GNSS Data Setting -----------------------%
% fpath_gnss: Path to GNSS data files. Do not include file names.
% fname_coord_gnss: file name of the site coordinates
% stn_id: station ID for sites to use for data-model comparison
% n_degree: n-th degree polynomial for detrending. 0 = remove mean, 1 = linear detrending, 2 = quadratic detrending
% use_berg_et_al: Logical, if true, GNSS rate will be adapted from Berg et al., 2024, if false, GNSS rate will be processed based on the raw data

use_berg_et_al = true;
fpath_gnss_new = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/gnss/DTU/';
fpath_gnss = '/Users/kyhan/Desktop/Projects/Greenland-LIA/data_obs/gnss/GNET_NEU_v2021_11_30/';
fname_coord_gnss = 'coordinates_v2021_11_30.txt';
%stn_1 = {'KAGZ','SCBY','KMOR','HRDG','JWLF','KMJP','JGBL','NORD','LEFN'};
%stn_2 = {'BLAS','NRSK','GROK','GMMA','YMER','DMHN','LBIB'};
%stn_3a = {'DANE','WTHG','HMBG','MSVG','DGJG','SCOR','VFDG'};
%stn_3b = {'KUAQ','MIK2'};
%stn_4 =  {'PLPK','KSNB','HEL2','KULU','ISOR','KBUG','LYNS','TREO','HJOR'};
%stn_5 = {'UTMG','TIMM','NNVN','QAQ1','SENU'}; % 'PAMI' has only 3 years of data. so exclude
%stn_6 = {'NUUK','KAPI','KLSQ','KELY','SISI','AASI','KAGA'};%,'ILUL','QEQE','QAAR','RINK'};
%stn_7 = {'SRMP','UPVK','KULL','ASKY','DKSG','THU2','MARG'};
%stn_id = [stn_1, stn_2, stn_3a, stn_3b, stn_4, stn_5, stn_6, stn_7];

stns = {'HRDG','JWLF','KMJP','JGBL','NORD','LEFN','BLAS','NRSK','GROK','GMMA','YMER','DMHN','LBIB','DANE','WTHG','HMBG','MSVG','DGJG','SCOR','VFDG','KUAQ','MIK2','PLPK','KSNB','HEL2','KULU','ISOR','KBUG','LYNS','TREO','HJOR','UTMG','TIMM','NNVN','QAQ1','SENU','NUUK','KAPI','KLSQ','KELY','SISI','AASI','KAGA','SRMP','UPVK','KULL','ASKY','DKSG','THU2','MARG','KAGZ','SCBY','KMOR'};
stn_id = stns;

n_degree = 0; % 0 = remove mean, 1 = linear detrending, 2 = quadratic detrending
annual_output = true; % preprocess data to be annual weighted mean
plot_gnss = true;
savefig_gnss = false;

%----------------------- GIA Data Setting -----------------------%
fpath_gia = '/Users/kyhan/Desktop/Projects/GIA-sensitivity-to-altimetry/data/gia/VLM_sites.dat';

%--------------------- Saltmarsh Data Setting --------------------%
% Setup for inputting saltmarsh data point coordinate
lat_ref = [63.470324, 62.547906, 63.131418]; %DM, TIMM, Vendom
lon_ref = [-41.926379, -42.27214, -41.457736];
fpath_rsl = '/Users/kyhan/Desktop/LIA-Project/Data/RSL_Woodroffe_2023/';
fname_rsl = 'SE_Greenland_SM_data.xlsx';