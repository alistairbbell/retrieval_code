%%%%%%%%% define forward model and retrieval settings %%%%%%%%%
%[Q, O] = setup_Q_O
%%% fixed rel!!!
function [Q, O] = setup_Q_O_ecmwf(Y,msm)

%=== Directories
top_dir = '/home/miawarac/retrievals/standard_retrievals/retrieval_data';

%arts_xmldata_path = atmlab( 'ARTS_XMLDATA_PATH' );
arts_includes     = atmlab( 'ARTS_INCLUDES' );
atmlab( 'VERBOSITY', 3);
atmlab( 'FMODEL_VERBOSITY', 0);

%%%%%%%%% initialize Q
Q=qarts;

%- General 
%
Q.INCLUDES            = { fullfile( arts_includes, 'general.arts'),fullfile( arts_includes, 'continua.arts' ) };
Q.ATMOSPHERE_DIM      = 1;
Q.STOKES_DIM          = 1;
Q.J_DO                = true;
Q.CLOUDBOX_DO         = false;
Q.INPUT_FILE_FORMAT   = 'double';
Q.OUTPUT_FILE_FORMAT  = 'ascii';

%- Radiative transfer
%
Q.Y_UNIT              = 'RJBT'; 
Q.YCALC_WSMS          = { 'yCalc' };
%
Q.PPATH_LMAX          = 250; % max distance between points
Q.PPATH_STEP_AGENDA   = { 'ppath_stepGeometric' };  % geom. Propagation, 
                                                    % refraction neglected

%- Surface
%
Q.R_GEOID            = constants( 'EARTH_RADIUS' );
Q.Z_SURFACE          = 16000;          % Just a dummy value. A 10 km
                                       % observation altitude is assumed here
Q.SENSOR_LOS         = Y.ZA;
Q.SENSOR_POS         = Q.Z_SURFACE;
Q.SENSOR_DO          = true;

%- Absorption
%
Q.ABS_LINES           = [top_dir,'/nedoluha1995_hyperfine.xml'];% fullfile( atmlab_example_data, 'o3line111ghz' );
Q.ABS_LINES_FORMAT    = 'Arts';
Q.ABS_LINESHAPE           = 'Voigt_Drayson';
Q.ABS_LINESHAPE_FACTOR    = 'quadratic';
Q.ABS_LINESHAPE_CUTOFF    = -1;
Q.ABSORPTION          = 'OnTheFly';
Q.ABS_NLS             = [];

%- Pressure grid
%
% z_toa                 = 110e3;
% Q.P_GRID              = z2p_simple( Q.Z_SURFACE-1e3 : 1000: z_toa )';
Q.P_GRID   = power(10,linspace(4.25,-2,100))';
Q.SENSOR_LOS          = Y.ZA;
Q.SENSOR_POS          = Y.Z_PLATFORM;
Q.SENSOR_DO           = true;


%- Frequency, spectrometer and pencil beam antenna
% The hypothetical spectrometer has rectangular response functions
%
Q.F_GRID              = Y.F;%read_datafile([top_dir,'/f_mono.aa'],'matrix');
%
H                     = qartsSensor;  % initialize sensor
%
H.SENSOR_NORM         = true;
%
% df                    = 0.5e6;
% H.F_BACKEND           = [ min(Q.F_GRID)+df : df : max(Q.F_GRID)-df ]';
H.F_BACKEND           = Y.F;
%
ch_resp = read_datafile(fullfile(top_dir,'channelresponse.aa'),'matrix');
B.name                = 'Spectrometer channel response function';
B.gridnames           = {'Frequency'} ;
B.grids               =  {ch_resp(:,1)'};
B.dataname            = 'Response';
B.data                = ch_resp(:,2)';
%
H.BACKEND_CHANNEL_RESPONSE{1} = B;
clear B ch_resp
%
Q.SENSOR_DO           = true;
Q.SENSOR_RESPONSE     = H;
%
Q.ANTENNA_DIM         = 1;
Q.MBLOCK_ZA_GRID      = 0;

% Antenna pattern    is included in correction for balancing/troposphere
% (a_trop)
% antenna_pattern = read_datafile(fullfile(curr_dir,'datafiles','miawara_antenna_pattern.aa'),'matrix');
% Q.SENSOR_RESPONSE.ANTENNA_DO = 1;
% Q.SENSOR_RESPONSE.ANTENNA_LOS = 0;
% Q.SENSOR_RESPONSE.ANTENNA_RESPONSE.name      = 'Antenna response';
% Q.SENSOR_RESPONSE.ANTENNA_RESPONSE.gridnames = {'Polarisation','Frequency','Zenith angle','Azimut angle'};%
% Q.SENSOR_RESPONSE.ANTENNA_RESPONSE.grids     = {{'1'},22.23508e9,antenna_pattern(:,1),0};
% Q.SENSOR_RESPONSE.ANTENNA_RESPONSE.dataname  = 'Response';
% Q.SENSOR_RESPONSE.ANTENNA_RESPONSE.data(1,1,:,1) = antenna_pattern(:,2);
% Q.MBLOCK_ZA_GRID  = antenna_pattern(:,1);
% clear antenna_pattern

%- Correlation of thermal noise (not correlated)
Q.TNOISE_C = covmat1d_from_cfun(Q.SENSOR_RESPONSE.F_BACKEND,[],'drc',0,0);  %only correlation, diagonal in Y.TNOISE

%%% PTZ

Q.HSE.ON       = false;
    mysql_open;
ptz = get_ptz_ECMWF_temporal_interpolation(datestr(msm.min_time,31), datestr(msm.max_time,31), msm.latitude, msm.longitude);
%ptz=get_pTz_ECMWF(msm.min_time,msm.max_time,msm.latitude,msm.longitude);
disp('Using ECMWF PTZ data')
%ptz=get_pT_AuraMLS_from_mysql_v2(msm.min_time,msm.max_time,msm.latitude,msm.longitude);

if isnan(nanmean(ptz(:,2)))
    Q = [];
    O = [];
    disp('No PTZ data found!')
    return
end

Q.T.RETRIEVE           = false;
Q.T.ATMDATA.TYPE       = 'atmdata';
Q.T.ATMDATA.NAME       = 'Temperature';
Q.T.ATMDATA.SOURCE     = 'MLS climatology';
Q.T.ATMDATA.DIM        = 1;
Q.T.ATMDATA.DATA       = ptz(:,2);
Q.T.ATMDATA.DATA_NAME  = 'Temperature';
Q.T.ATMDATA.DATA_UNIT  = 'K';
Q.T.ATMDATA.GRID1      = ptz(:,1);
Q.T.ATMDATA.GRID1_NAME = 'Pressure';
Q.T.ATMDATA.GRID1_UNIT = 'Pa';
Q.Z_ATMDATA.TYPE       = 'atmdata';
Q.Z_ATMDATA.NAME       = 'Altitude';
Q.Z_ATMDATA.SOURCE     = 'MLS';
Q.Z_ATMDATA.DIM        = 1;
Q.Z_ATMDATA.DATA       = ptz(:,3);
Q.Z_ATMDATA.DATA_NAME  = 'Altitude';
Q.Z_ATMDATA.DATA_UNIT  = 'm';
Q.Z_ATMDATA.GRID1      = ptz(:,1);
Q.Z_ATMDATA.GRID1_NAME = 'Pressure';
Q.Z_ATMDATA.GRID1_UNIT = 'Pa';
clear ptz

%- Species
% Water vapor
Q.ABS_SPECIES(1).TAG      = {'H2O'};
Q.ABS_SPECIES(1).RETRIEVE = true;
Q.ABS_SPECIES(1).L2       = true;
Q.ABS_SPECIES(1).GRIDS{1} = power(10,linspace(4,-2,50))'; %somewhat denser than old grid
Q.ABS_SPECIES(1).GRIDS{2} = []; %lat grid
Q.ABS_SPECIES(1).GRIDS{3} = []; %lon grid
Q.ABS_SPECIES(1).UNIT     = 'vmr';

%%% a priori profile
h2o=get_apriori_h2o_from_mysql(datenum(Y.YEAR, Y.MONTH, Y.DAY, Y.HOUR, Y.MINUTE, Y.SECOND),Y.LATITUDE);

Q.ABS_SPECIES(1).ATMDATA.TYPE       = 'atmdata';
Q.ABS_SPECIES(1).ATMDATA.NAME       = {'H2O'};
Q.ABS_SPECIES(1).ATMDATA.SOURCE     = 'mls climatology';
Q.ABS_SPECIES(1).ATMDATA.DIM        = 1;
Q.ABS_SPECIES(1).ATMDATA.DATA       = h2o(:,2);
Q.ABS_SPECIES(1).ATMDATA.DATA_NAME  = '';
Q.ABS_SPECIES(1).ATMDATA.DATA_UNIT  = '-';      %%% ???rel
Q.ABS_SPECIES(1).ATMDATA.GRID1      = h2o(:,1);
Q.ABS_SPECIES(1).ATMDATA.GRID1_NAME = 'Pressure';
Q.ABS_SPECIES(1).ATMDATA.GRID1_UNIT = 'Pa';


% read apriori covariance matrix relative and absolute
sx_ppm=read_datafile([top_dir,'/sx.apriori.H2O.ppm.aa'],'matrix');

% convert sigma from ppm to fraction of a priori profile, on a priori grid
% sx_rel(:,1)=Q.ABS_SPECIES(1).ATMDATA.GRID1;
% sx_rel(:,2)=interp1(log10(sx_ppm(:,1)),sx_ppm(:,2),log10(Q.ABS_SPECIES(1).ATMDATA.GRID1))./h2o(:,2);
% 
% % constrain the relative sigma to 15-65%
% sx_rel(sx_rel(:,2)<0.15 | sx_rel(:,1)>4,2)=0.15;
% sx_rel(sx_rel(:,2)>0.65,2)=0.65;
% 
% % smooth the relative sigma profile
% sx_rel(:,2)=smooth(sx_rel(:,2),10);
% 
% %%sx_rel=read_datafile([top_dir,'/sx.apriori.H2O.rel.aa'],'matrix');    sx
% %%with fixed rel
% 
% sx_abs=[sx_rel(:,1),sx_rel(:,2).*Q.ABS_SPECIES(1).ATMDATA.DATA];        % must match Q.ABS_SPECIES.UNIT
sx_abs(:,1)=Q.ABS_SPECIES(1).ATMDATA.GRID1;
sx_abs(:,2)=interp1(log10(sx_ppm(:,1)),sx_ppm(:,2),log10(Q.ABS_SPECIES(1).ATMDATA.GRID1));
corr_length=0.25;

% write it to file
Q.ABS_SPECIES(1).SX = covmat1d_from_cfun( Q.ABS_SPECIES(1).GRIDS{1}, sx_abs, 'exp', corr_length, 0, @log10);


%- Polyfit
%
% A polynomial of order 6 is used for "baseline fit".
%
Q.POLYFIT.RETRIEVE        = true;
Q.POLYFIT.ORDER           = 6;
Q.POLYFIT.L2              = true;
Q.POLYFIT.SX0             = 1^2; 
Q.POLYFIT.SX1             = 0.5^2; 
Q.POLYFIT.SX2             = 0.2^2;
Q.POLYFIT.SX3             = 0.1^2;
Q.POLYFIT.SX4             = 0.1^2; 
Q.POLYFIT.SX5             = 0.1^2; 
Q.POLYFIT.SX6             = 0.1^2; 

%- Define L2 structure (beside retrieval quantities below)
%  Q.L2_EXTRA = {'dx','cost','e','eo','es','yf','A','S','So','Ss','G','J','mresp','date','xa','y','bl','tnoise','ptz'};
Q.L2_EXTRA = {'dx','cost','e','eo','es','yf','A','S','So','Ss','mresp','date','xa','y','bl','tnoise','ptz'};




% %- Frequency (shift retrieval here, shift+stretch also possible)
% %
% Q.FFIT.RETRIEVE           = true;
% Q.FFIT.DF                 = 50e3; 
% Q.FFIT.ORDER              = 0;
% Q.FFIT.SX                 = 50e3^2;  
% Q.FFIT.L2                 = true;
% 
% Q.FFIT.RETRIEVE           = true;
% Q.FFIT.DF                 = 100e3; 
% Q.FFIT.ORDER              = 1;
% Q.FFIT.SX                 = [300e3^2 0;0 100e3^2]; 
% Q.FFIT.L2                 = true;

%%% Define O
O = qp2_l2( Q );  % This ensures that OEM returns the varibles needed
O.linear = true;
