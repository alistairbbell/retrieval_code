% function Q = createQ_with_sideband(Y,R)
%
% This function creates the Qpack2 Q-structure, based on
% information about the measurements (Y) and retrieval specific
% information (R).
%
%
% DS, 2010-06-30
%
function Q = createQ_with_sideband(Y,R)

% Initialize Atmlab settings
arts_xmldata_path = atmlab( 'ARTS_XMLDATA_PATH' );
arts_includes     = atmlab( 'ARTS_INCLUDES' );
atmlab( 'VERBOSITY', 3);
atmlab( 'FMODEL_VERBOSITY', 0);
fascod = fullfile( arts_xmldata_path, 'atmosphere', 'fascod' );
curr_dir = pwd;

% Initilize the Q-Structure
Q = qarts;

%- General
Q.INCLUDES            = { fullfile( arts_includes, 'general.arts' ), fullfile( arts_includes, 'continua.arts' ), fullfile( curr_dir,'datafiles', 'continua_h2o.arts')};
Q.INPUT_FILE_FORMAT   = 'double';
Q.OUTPUT_FILE_FORMAT  = 'ascii';
Q.ATMOSPHERE_DIM      = 1;
% Q.USE_RAW_ATMOSPHERE  = false;
Q.HSE.ON              = false;
Q.CLOUDBOX_DO         = false;
Q.STOKES_DIM          = 1;
Q.J_DO                = true;

%- Radiative transfer
Q.Y_UNIT              = 'RJBT';
% Q.WSMS_BEFORE_RTE     = { 'atm_checkedCalc' };
Q.YCALC_WSMS          = { 'yCalc' };
Q.PPATH_LMAX          = 200;
Q.PPATH_STEP_AGENDA   = { 'ppath_stepGeometric' };

%- Absorption
Q.ABS_LINES            = fullfile(curr_dir,'datafiles','lines_hitran06_jpl01_08.xml');
Q.ABS_LINES_FORMAT     = 'Arts';
Q.ABS_LINESHAPE        = 'Voigt_Kuntz6'; % See 'arts2 -d abs_lineshapeDefine' for more line shapes
Q.ABS_LINESHAPE_CUTOFF = -1;
Q.ABS_LINESHAPE_FACTOR = 'VVH';
Q.ABSORPTION           = 'OnTheFly';
Q.ABS_NLS              = []; % See 'arts2 -d abs_nls'


%- Frequency, spectrometer and antenna
Q.F_GRID              = read_datafile(fullfile(curr_dir,'datafiles','miawara_opt_fmono_ext.aa'),'matrix');
Q.SENSOR_LOS          = Y.ZA;
Q.SENSOR_POS          = Y.Z_PLATFORM;
Q.SENSOR_DO           = true;

% Initialite the SENSOR_RESPONSE structure
Q.SENSOR_RESPONSE     = qartsSensor;
Q.SENSOR_RESPONSE.SENSOR_NORM  = true;

% Acqiris FFT channel response
ch_resp = read_datafile(fullfile(curr_dir,'datafiles','FFT_ACQ1_channelresponse.aa'),'matrix');
Q.SENSOR_RESPONSE.BACKEND_DO = 1;
Q.SENSOR_RESPONSE.F_BACKEND  = Y.F;
Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.name      = 'Acqiris FFT channel response';
Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.gridnames = {'Frequency'};
Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.grids     = {ch_resp(:,1)'};
Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.dataname  = 'Response';
Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.data      = ch_resp(:,2)';
clear ch_resp

Q.ANTENNA_DIM     = 1;
Q.MBLOCK_AA_GRID  = {};
Q.MBLOCK_ZA_GRID  = 0;

% Antenna pattern      Remark: Antenna pattern should not be considered, because with balancing, the normal antenna pattern does not apply anymore!
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

% Sideband response
Q.SENSOR_RESPONSE.LO = 20.135e9;
Q.SENSOR_RESPONSE.F_BACKEND = Q.SENSOR_RESPONSE.F_BACKEND - Q.SENSOR_RESPONSE.LO;
sb_response = read_datafile(fullfile(curr_dir,'datafiles','sbfilter.aa'),'matrix');
sb_response_ip(:,1) = linspace(min(Q.F_GRID),max(Q.F_GRID),100);
sb_response_ip(:,2) = interp1(sb_response(:,1),sb_response(:,2),sb_response_ip(:,1));
sb_response_ip(:,1) = sb_response_ip(:,1)-Q.SENSOR_RESPONSE.LO;
Q.SENSOR_RESPONSE.MIXER_DO = 1;
Q.SENSOR_RESPONSE.SIDEBAND_RESPONSE.name      = 'Sideband and Mixer response';
Q.SENSOR_RESPONSE.SIDEBAND_RESPONSE.gridnames = {'Frequency'};
Q.SENSOR_RESPONSE.SIDEBAND_RESPONSE.grids     = {sb_response_ip(:,1)'};
Q.SENSOR_RESPONSE.SIDEBAND_RESPONSE.dataname  = 'Response';
Q.SENSOR_RESPONSE.SIDEBAND_RESPONSE.data      = sb_response_ip(:,2)';
clear sb_response sb_response_ip

%- Correlation of thermal noise (not correlated)
Q.TNOISE_C = covmat1d_from_cfun(Q.SENSOR_RESPONSE.F_BACKEND,[],'drc',0,0);

% Temperature and altitude data
% Q.T.ATMDATA = gf_artsxml( fullfile(arts_xmldata_path,'atmosphere','fascod','midlatitude-winter.t.xml'),'Temperature','t_field');
% Q.Z_ATMDATA = gf_artsxml( fullfile(arts_xmldata_path,'atmosphere','fascod','midlatitude-winter.z.xml'),'Altitude','z_field');

if isfield(R,'PTZ_profile_hour_interval') && ~isempty(R.PTZ_profile_hour_interval) && ~strcmp(R.PTZ_profile,'ECMWF')
    error('If hour_interval is set, PTZ profile must be ECMWF!')
elseif ~isfield(R,'PTZ_profile_hour_interval')
    R.PTZ_profile_hour_interval = [];
end
switch R.PTZ_profile
    case 'ECMWF'
        ptz = get_ptz_ECMWF(datestr(Y.MIN_TIME,31),datestr(Y.MAX_TIME,31),Y.LATITUDE,Y.LONGITUDE,R.PTZ_profile_hour_interval);
    case 'USSTD'
        ptz_raw = [[101325 22632.1 5474.89 868.019 110.906 66.9389 3.95642 0.3734 0.01]',[288.15 216.65 216.65 228.65 270.65 270.65 214.65 273.15-86.2 273.15-86.2]',[0 11000 20000 32000 47000 51000 71000 84852 100000]'];
        ptz(:,1) = power(10,linspace(log10(ptz_raw(1,1)),log10(ptz_raw(end,1)),100))';
        ptz(:,2) = interp1(log10(ptz_raw(:,1)),ptz_raw(:,2),log10(ptz(:,1)));
        ptz(:,3) = interp1(log10(ptz_raw(:,1)),ptz_raw(:,3),log10(ptz(:,1)));
        ptz(:,2) = smooth(ptz(:,2),5);
        ptz(:,3) = smooth(ptz(:,3),5);
    case 'MLS_actual'
        ptz = get_ptz_MLS_from_mysql(datestr(Y.MIN_TIME,31),datestr(Y.MAX_TIME,31),Y.LATITUDE,Y.LONGITUDE,3.3);
    case 'MLS_climat'
        ptz = get_ptz_MLSclimat(mean([Y.MIN_TIME,Y.MAX_TIME]));
    case 'MLS_climat_2004_2010'
        ptz = get_ptz_MLSclimat_2004_2010(mean([Y.MIN_TIME,Y.MAX_TIME]));
    otherwise
        error('R.PTZ_profile was not understood!')
end
if isempty(ptz)
    Q = [];
    disp('No PTZ data found!')
    return
end
Q.T.RETRIEVE           = false;
Q.T.ATMDATA.TYPE       = 'atmdata';
Q.T.ATMDATA.NAME       = 'Temperature';
Q.T.ATMDATA.SOURCE     = R.PTZ_profile;
Q.T.ATMDATA.DIM        = 1;
Q.T.ATMDATA.DATA       = ptz(:,2);
Q.T.ATMDATA.DATA_NAME  = 'Temperature';
Q.T.ATMDATA.DATA_UNIT  = 'K';
Q.T.ATMDATA.GRID1      = ptz(:,1);
Q.T.ATMDATA.GRID1_NAME = 'Pressure';
Q.T.ATMDATA.GRID1_UNIT = 'Pa';
Q.Z_ATMDATA.TYPE       = 'atmdata';
Q.Z_ATMDATA.NAME       = 'Altitude';
Q.Z_ATMDATA.SOURCE     = R.PTZ_profile;
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
Q.ABS_SPECIES(1).TAG      = { R.H2O_model };
Q.ABS_SPECIES(1).RETRIEVE = true;
Q.ABS_SPECIES(1).L2       = true;
Q.ABS_SPECIES(1).GRIDS{1} = power(10,linspace(4,-2,50))';
Q.ABS_SPECIES(1).GRIDS{2} = [];
Q.ABS_SPECIES(1).GRIDS{3} = [];
Q.ABS_SPECIES(1).UNIT     = 'vmr';
switch R.H2O_apriori_profile
    case 'USSTD'
        h2odata = read_datafile(fullfile(curr_dir,'datafiles','H2O_USSTD.aa'),'matrix');
    case 'MLS_mean'
        
    case 'MLS_climat'
        h2odata = get_h2o_apr_MLSclimat(datenum(Y.YEAR, Y.MONTH, Y.DAY, Y.HOUR, Y.MINUTE, Y.SECOND),Y.LATITUDE);
    case 'custom'
        h2odata = R.H2O_apriori_profile_data;
    otherwise
        error('R.H2O_apriori_profile was not understood!')
end
Q.ABS_SPECIES(1).ATMDATA.TYPE       = 'atmdata';
Q.ABS_SPECIES(1).ATMDATA.NAME       = 'H2O';
Q.ABS_SPECIES(1).ATMDATA.SOURCE     = R.H2O_apriori_profile;
Q.ABS_SPECIES(1).ATMDATA.DIM        = 1;
Q.ABS_SPECIES(1).ATMDATA.DATA       = h2odata(:,2);
Q.ABS_SPECIES(1).ATMDATA.DATA_NAME  = 'Volume mixing ratio';
Q.ABS_SPECIES(1).ATMDATA.DATA_UNIT  = '-';
Q.ABS_SPECIES(1).ATMDATA.GRID1      = h2odata(:,1);
Q.ABS_SPECIES(1).ATMDATA.GRID1_NAME = 'Pressure';
Q.ABS_SPECIES(1).ATMDATA.GRID1_UNIT = 'Pa';
switch R.H2O_apriori_covariance
    case 'v17'
        SX_rel = read_datafile(fullfile(curr_dir,'datafiles','v17.sx.apriori.H2O.aa'),'aomatrix');
    case 'v18'
        SX_rel = read_datafile(fullfile(curr_dir,'datafiles','v18.sx.apriori.H2O.aa'),'aomatrix');
    case 'v22'
        SX_v22_abs = read_datafile(fullfile(curr_dir,'datafiles','v22.sx.apriori.H2O.ppm.aa'),'matrix');
        SX_abs{1} = [3;0];
        SX_abs{2}(:,1) = log10(Q.ABS_SPECIES(1).ATMDATA.GRID1);
        SX_abs{2}(:,2) = interp1(SX_v22_abs(:,1),SX_v22_abs(:,2),Q.ABS_SPECIES(1).ATMDATA.GRID1);
        SX_abs{2}(:,3) = 0.25;
        rel = SX_abs{2}(:,2)./Q.ABS_SPECIES(1).ATMDATA.DATA;
        rel(rel<0.15 | rel>4) = 0.15;
        rel(rel>0.65) = 0.65;
        rel(isnan(rel)) = 0.65;
        rel(1) = rel(2);
        rel=smooth(rel,10);
        SX_abs{2}(:,2) = rel.*Q.ABS_SPECIES(1).ATMDATA.DATA;
        clear rel SX_v22_abs
    otherwise
        error('R.H2O_apriori_covariance was not understood!')
end
if exist('SX_rel','var')
    SX_abs{1} = SX_rel{1};
    SX_abs{2}(:,1) = log10(Q.ABS_SPECIES(1).ATMDATA.GRID1);
    SX_abs{2}(:,2) = interp1(SX_rel{2}(:,1),SX_rel{2}(:,2),SX_abs{2}(:,1)).*Q.ABS_SPECIES(1).ATMDATA.DATA;
    SX_abs{2}(:,3) = interp1(SX_rel{2}(:,1),SX_rel{2}(:,3),SX_abs{2}(:,1));
elseif ~exist('SX_abs','var')
    error('No H2O covariance matrix defined!')
end
switch SX_abs{1}(1)
    case 0
        corr_fun = 'drc';
    case 1
        corr_fun = 'linn';
    case 2
        corr_fun = 'exp';
    case 3
        corr_fun = 'gau';
    otherwise
        error('Correlation function of the H2O covariance matrix not understood!')
end
Q.ABS_SPECIES(1).SX = covmat1d_from_cfun( Q.ABS_SPECIES(1).GRIDS{1}, [power(10,SX_abs{2}(:,1)) SX_abs{2}(:,2)], corr_fun, [power(10,SX_abs{2}(:,1)) SX_abs{2}(:,3)], SX_abs{1}(2), @log10);
clear SX_rel SX_abs h2odata


%- Pressure grid (not the retrieval grid!)
% Q.P_GRID   = power(10,linspace(max(log10(Q.ABS_SPECIES(1).GRIDS{1})), min(log10(Q.ABS_SPECIES(1).GRIDS{1})) , 2*length(Q.ABS_SPECIES(1).GRIDS{1})))';
Q.P_GRID   = power(10,linspace(4.25,-2,100))';

%- Surface
Q.R_GEOID            = constants( 'EARTH_RADIUS' );
Q.Z_SURFACE          = 16000; %interp1(log10(Q.Z_ATMDATA.GRID1),Q.Z_ATMDATA.DATA,max(log10(Q.P_GRID)));
Q.SENSOR_LOS         = Y.ZA;
Q.SENSOR_POS         = Q.Z_SURFACE;
Q.SENSOR_DO          = true;



%- Define L2 structure (beside retrieval quantities below)
Q.L2_EXTRA = {'dx','cost','e','eo','es','yf','A','S','So','Ss','mresp','date','xa','y','bl','tnoise','ptz'};

