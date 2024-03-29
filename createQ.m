%Q = createQ(Y,R, ptz)

% The Qpack2.4 structure is created in here. 
% 
% Here ptz is contains the pressure to altitude relation where the first
% column is the pressure (in Pa), the second column is temperature (in K), 
% and the third column is the altitude (in m).
% 
% 
% To note that some fields will not be compatible
% with precvious versions of ARTS and so ARTS 2.4 and ATMLAB 
% 2.4 should be used if using this routine.

% 2022-08-20
%
function Q = createQ(Y,R, ptz)

% Initialize Atmlab settings
arts_xmldata_path = atmlab( 'ARTS_XMLDATA_PATH' );
arts_includes     = atmlab( 'ARTS_INCLUDES' );
atmlab( 'VERBOSITY', 3);
atmlab( 'FMODEL_VERBOSITY', 0);
curr_dir ='/home/alistair/MIAWARA_ret/old_miawara_calibration/retrievals/retrievals_old/';
disp('curr_dir')
disp(curr_dir)

% Initilize the Q-Structure
Q = qarts;

%- General
Q.INCLUDES            = { fullfile( arts_includes, 'general.arts' ), fullfile( arts_includes, 'continua.arts' ), fullfile( curr_dir,'datafiles', 'continua_h2o.arts')};
Q.INPUT_FILE_FORMAT   = 'double';
Q.OUTPUT_FILE_FORMAT  = 'ascii';
Q.ATMOSPHERE_DIM      = 1; % Atmospheric dimensionality (1-3)
Q.HSE.ON              = false; % Boolean to activate the hydrostatic equilibrium
Q.CLOUDBOX_DO         = false; % Boolean to activate the cloud box
Q.STOKES_DIM          = 1; % Dimensionality of the Stokes vector (1-4)
Q.J_DO                = true; % Boolean to include calculation of Jacobians

%- Radiative Transfer
%Q.Y_UNIT              = 'RJBT'; % Radiance units (1, RJBT or PlanckBT)
Q.YCALC_WSMS          = { 'yCalc' };
Q.PPATH_LMAX          = 200;
Q.PPATH_STEP_AGENDA   = { 'ppath_stepGeometric' };

%- Absorption
%Q.ABS_LINES            = fullfile(curr_dir,'datafiles','lines_hitran06_jpl01_08_finesplitting.xml');
%Q.ABS_LINES            = fullfile(curr_dir,'datafiles','lines_hitran06_jpl01_08_finesplitting_only_H2O.xml');
%Q.ABS_LINES            = fullfile(curr_dir,'datafiles','lines_nedoluha1995_finesplitting.xml');
%AB Q.ABS_LINES            = fullfile(curr_dir,'datafiles','lines_nedoluha1995_finesplitting_only_H2O.xml');
%AB Q.ABS_LINES_FORMAT     = 'Arts';
%AB Q.ABS_LINESHAPE        = 'Voigt_Kuntz6'; % See 'arts2 -d abs_lineshapeDefine' for more line shapes

%disp('lineshapes')
%arts -d abs_lineshapeDefine

%AB Q.ABS_LINESHAPE_CUTOFF = -1;
%AB Q.ABS_LINESHAPE_FACTOR = 'VVH';
Q.ABSORPTION           = 'OnTheFly'; %on the fly must be used unless lookup tables are created in advance
Q.ABS_NLS              = []; % See 'arts2 -d abs_nls'


%- Frequency, spectrometer and antenna
%Q.F_GRID              = read_datafile(fullfile(curr_dir,'datafiles','miawara_opt_fmono_ext.aa'),'matrix');

Q.SENSOR_LOS          = Y.ZA;
Q.SENSOR_POS          = Y.Z_PLATFORM;
Q.SENSOR_DO           = false; 
% sensor_do - Boolean to include sensor characteristics. Otherwise monochromatic
%pencil beam calculations are performed. (los and sensor pos always used)

% Initialite the SENSOR_RESPONSE structure
 %Q.SENSOR_RESPONSE     = qartsSensor;
%Q.SENSOR_RESPONSE.SENSOR_NORM  = true;

% Acqiris FFT channel response
ch_resp = read_datafile(fullfile(curr_dir,'datafiles','FFT_ACQ1_channelresponse.aa'),'matrix');
Q.SENSOR_RESPONSE.BACKEND_DO = 1;
Q.SENSOR_RESPONSE.MIXER_DO = 0;
Q.SENSOR_RESPONSE.F_BACKEND  = Y.F;
Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.name      = 'Acqiris FFT channel response';
Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.gridnames = {'Frequency'};
Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.grids     = {ch_resp(:,1)'};
Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.dataname  = 'Response';
Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.data      = ch_resp(:,2)';
clear ch_resp

Q.ANTENNA_DIM     = 1;
%Q.MBLOCK_AA_GRID  = {};
%Q.MBLOCK_ZA_GRID  = 0;

%- Correlation of thermal noise (not correlated)
Q.TNOISE_C = covmat1d_from_cfun(Q.SENSOR_RESPONSE.F_BACKEND,[],'drc',0,0);

% Temperature and altitude data
% Q.T.ATMDATA = gf_artsxml( fullfile(arts_xmldata_path,'atmosphere','fascod','midlatitude-winter.t.xml'),'Temperature','t_field');
% Q.Z_ATMDATA = gf_artsxml( fullfile(arts_xmldata_path,'atmosphere','fascod','midlatitude-winter.z.xml'),'Altitude','z_field');
if isempty(ptz)
    Q = [];
    disp('No PTZ data found!')
    return
end

clear ptz


Q.SENSOR_LOS         = Y.ZA;
Q.SENSOR_POS         = Q.Z_SURFACE;
Q.SENSOR_DO          = true;

%- Define L2 structure (beside retrieval quantities below)
Q.L2_EXTRA = {'dx','cost','e','eo','es','yf','A','S','So','Ss','mresp','date','xa','y','bl','tnoise','ptz'};

