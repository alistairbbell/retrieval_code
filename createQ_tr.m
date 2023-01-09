% function Q = createQ(Y,R)
%
% This function creates the Qpack2 Q-structure, based on
% information about the measurements (Y) and retrieval specific
% information (R).
%
%
% DS, 2010-06-30
% Adapted for trop. retrieval, RB, 2010-07-19
%
function [Q,R] = createQ_tr(Y,R)
mtime    = datenum([Y.YEAR Y.MONTH Y.DAY Y.HOUR Y.MINUTE Y.SECOND]);
curr_dir = '/home/miawara/retrievals/retrieval_arts2_troposphere/';

% Create the Q-Structure required by Qpack2 and set some basic variables
Q = qarts;

%- General
arts_includes          = atmlab( 'ARTS_INCLUDES' );
atmlab( 'VERBOSITY', 3);
atmlab( 'FMODEL_VERBOSITY', 0);

Q.INCLUDES             = { fullfile( arts_includes, 'general.arts' ), fullfile( arts_includes, 'continua.arts' )};
Q.INPUT_FILE_FORMAT    = 'double';
Q.OUTPUT_FILE_FORMAT   = 'ascii';
Q.R_GEOID              = constants( 'EARTH_RADIUS' );
Q.Z_SURFACE            = Y.Z_PLATFORM;
Q.ATMOSPHERE_DIM       = 1;        % Atmospheric dimensionality (1-3)
Q.HSE.ON               = false;    % Boolean to activate the hydrostatic equilibrium
Q.CLOUDBOX_DO          = false;    % Boolean to activate the cloud box
Q.STOKES_DIM           = 1;        % Dimensionality of the Stokes vector (1-4)
Q.J_DO                 = true;     % Boolean to include calculation of Jacobians

%- Radiative transfer
Q.Y_UNIT               = 'RJBT';   % Radiance units (1, RJBT or PlanckBT)
Q.WSMS_BEFORE_RTE      = { 'basics_checkedCalc' };
Q.YCALC_WSMS           = { 'yCalc' };
Q.PPATH_LMAX           = 200;
Q.PPATH_STEP_AGENDA    = { 'ppath_stepGeometric' };

%- Absorption
Q.ABS_LINES            = fullfile(curr_dir,'datafiles','lines_hitran06_jpl01_08.xml');
Q.ABS_LINES_FORMAT     = 'Arts';
Q.ABS_LINESHAPE        = 'Lorentz'; % ['Voigt_Kuntz6'], See 'arts2 -d abs_lineshapeDefine' for more line shapes
Q.ABS_LINESHAPE_CUTOFF = -1;
Q.ABS_LINESHAPE_FACTOR = 'quadratic'; % ['VVH']
Q.ABSORPTION           = 'OnTheFly';
Q.ABS_NLS              = []; % See 'arts2 -d abs_nls'


%- Frequency, spectrometer and antenna
Q.F_GRID              = xmlLoad(fullfile(curr_dir,'datafiles','f_grid_sb.xml'));
Q.SENSOR_LOS          = Y.ZA;
Q.SENSOR_POS          = Y.Z_PLATFORM;
Q.SENSOR_DO           = true;

% Initialize the SENSOR_RESPONSE structure
Q.SENSOR_RESPONSE     = qartsSensor;
Q.SENSOR_RESPONSE.SENSOR_NORM  = true;

% Acqiris FFT channel response
if isfield(R,'ch_resp')
    ch_resp = R.ch_resp;
else
    ch_resp               = xmlLoad(fullfile(curr_dir,'datafiles','channelresponse_20MHz.xml'));
end

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

% Sideband response
Q.SENSOR_RESPONSE.LO = 20.135e9;
Q.SENSOR_RESPONSE.F_BACKEND = Q.SENSOR_RESPONSE.F_BACKEND - Q.SENSOR_RESPONSE.LO;
sb_response = xmlLoad(fullfile(curr_dir,'datafiles','sbfilter.xml'));
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

%- Correlation of thermal noise
Q.TNOISE_C = covmat1d_from_cfun(Q.SENSOR_RESPONSE.F_BACKEND,[],'exp',6E5,0.3); % ['drc',0,0]

% Temperature and altitude data
switch R.PTZ_profile
    case 'ECMWF'
        ptz = get_ptz_ECMWF_temporal_interpolation(datestr(Y.MIN_TIME,31),datestr(Y.MAX_TIME,31),Y.LATITUDE,Y.LONGITUDE);
    case 'ECMWF_nodiurnal'  
        ptz = get_ptz_ECMWF_temporal_interpolation(datestr(mean([Y.MIN_TIME Y.MAX_TIME])-1,31),datestr(mean([Y.MIN_TIME Y.MAX_TIME])+1,31),Y.LATITUDE,Y.LONGITUDE);
    case 'ECMWF_trop'
        [ptz,vmr] = get_pTz_ECMWF_trop(mtime-0.5,mtime+0.5,47,7,1);                  
    case 'ECMWF+MLS'
        ptz_mls   = get_ptz_MLS_from_mysql(datestr(mtime-0.5,31),datestr(mtime+0.5,31),Y.LATITUDE,Y.LONGITUDE,3.3);
        ptz_ecmwf = get_pTz_ECMWF_trop(mtime-0.5,mtime+0.5,47,7,1); 
        
        if isempty(ptz_ecmwf)||isempty(ptz_mls)
            Q = [];
            return       
        end
        ptz       = [ptz_ecmwf(ptz_ecmwf(:,1)>10e2,:); ptz_mls(ptz_mls(:,1)<10e2,:)];
    case 'USSTD'
        ptz_raw = [[101325 22632.1 5474.89 868.019 110.906 66.9389 3.95642 0.3734 0.01]',[288.15 216.65 216.65 228.65 270.65 270.65 214.65 273.15-86.2 273.15-86.2]',[0 11000 20000 32000 47000 51000 71000 84852 100000]'];
        ptz(:,1) = power(10,linspace(log10(ptz_raw(1,1)),log10(ptz_raw(end,1)),100))';
        ptz(:,2) = interp1(log10(ptz_raw(:,1)),ptz_raw(:,2),log10(ptz(:,1)));
        ptz(:,3) = interp1(log10(ptz_raw(:,1)),ptz_raw(:,3),log10(ptz(:,1)));
        ptz(:,2) = smooth(ptz(:,2),5);
        ptz(:,3) = smooth(ptz(:,3),5);
    case 'MLSv3.3'
        ptz = get_ptz_MLS_from_mysql(datestr(Y.MIN_TIME,31),datestr(Y.MAX_TIME,31),Y.LATITUDE,Y.LONGITUDE,3.3);
    case 'MLSv2.2'
        ptz = get_ptz_MLS_from_mysql(datestr(Y.MIN_TIME,31),datestr(Y.MAX_TIME,31),Y.LATITUDE,Y.LONGITUDE,2.2);
    case 'MLS_climatv2.2'
        ptz = get_ptz_MLSclimat(mean([Y.MIN_TIME,Y.MAX_TIME]));
    case 'MLS_climat_v3.3'
        ptz = get_ptz_MLSclimat_2004_2010(mean([Y.MIN_TIME,Y.MAX_TIME]));
    otherwise
        error('R.PTZ_profile was not understood!')
end
if isempty(ptz)
    Q = [];
    disp('No ptz data found!')
    return
end

if isfield(R,'surface')
    switch R.surface            
        case 'meteo'
            if (Y.MAX_TIME-Y.MIN_TIME)<10/60/24
                met = get_miawara_meteo(datestr(Y.MIN_TIME-5/60/24,31),...
                    datestr(Y.MAX_TIME+5/60/24,31));
            else
                met = get_miawara_meteo(datestr(Y.MIN_TIME,31),...
                    datestr(Y.MAX_TIME,31));
            end
            if isempty(met)
                Q = [];
                disp('No meteo data found');
                return
            end
            ptz_ground = [nanmean(met(:,2)) nanmean(met(:,3)) 905];
            vmr_ground = rh2vmr(ptz_ground(1),ptz_ground(2),nanmean(met(:,4)));
            
            R.surf_data = [ptz_ground vmr_ground nanmean(met(:,4))];
        case 'snd'
            ptz_ground = R.ptz_ground;
            vmr_ground = R.vmr_ground;
        otherwise
            error('Surface value type not known');
    end
    
    % Merge surface data and ptz-data
    ind = find(ptz(:,3)>905 & ptz(:,1)<ptz_ground(1));
    ptz = [ptz_ground;ptz(ind,:)];
    if exist('vmr','var')
        vmr = [vmr_ground;vmr(ind)];
    else
        vmr = vmr_ground;
    end
end
if strcmp(R.PTZ_profile,'ECMWF_trop')||strcmp(R.PTZ_profile,'ECMWF+MLS')
    R.ptz_ecmwf = ptz;
    R.vmr_ecmwf = vmr;
end
            
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
Q.T.RETRIEVE           = false;
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

%- Pressure grid (not the retrieval grid!)
%Q.P_GRID               = power(10,linspace(5,-2,100))';
% p-grid (standard grid + surface pressure)
if isfield(R,'P_GRID_calc')
    Q.P_GRID               = [ptz(1,1);R.P_GRID_calc(R.P_GRID_calc<ptz(1,1))];
else
    p1                     = logspace(5,4,17);
    p2                     = logspace(4,1,13);
    p                      = [p1 p2(2:end)]';
    Q.P_GRID               = [ptz(1,1);p(p<ptz(1,1))];
end


%- Species
% Water vapor
Q.ABS_SPECIES(1).TAG      = { R.H2O_model };
Q.ABS_SPECIES(1).RETRIEVE = true;
Q.ABS_SPECIES(1).L2       = true;
if isfield(R,'P_GRID_ret')
    Q.ABS_SPECIES(1).GRIDS{1} = [ptz(1,1);R.P_GRID_ret(R.P_GRID_ret<ptz(1,1))];
else
    Q.ABS_SPECIES(1).GRIDS{1} = Q.P_GRID;
end
Q.ABS_SPECIES(1).GRIDS{2} = [];
Q.ABS_SPECIES(1).GRIDS{3} = [];
Q.ABS_SPECIES(1).UNIT     = R.unit;

switch R.H2O_apriori_profile
    case 'USSTD'
        h2odata = xmlLoad(fullfile(curr_dir,'datafiles','H2O_USSTD.xml'));        
    case 'MLS_climat'
        h2odata = get_h2o_apr_MLSclimat(mtime,Y.LATITUDE);
    case 'RS_THUN'     
        h2odata = get_apriori_h2o_from_RS_THUN(datenum([Y.YEAR Y.MONTH Y.DAY Y.HOUR Y.MINUTE Y.SECOND]),[ptz(1,1) vmr(1)]); % RS-THUN + surface value (Zimm-meteo)
                
    case 'RS_THUN-MIA'
        xa  = get_apriori_h2o_from_RS_THUN(mtime,[ptz(1,1) vmr(1)]); % RS-THUN + surface value (Zimm-meteo)
        mia = get_miawara_level2(datestr(floor(mtime),'yyyy-mm-dd'),datestr(ceil(mtime),'yyyy-mm-dd'),'24');
        
        if ~isempty(mia)
            h2odata = [xa(xa(:,1)>17000,:); [mia.pressure mia.vmr]];
            clear indt mia
        else
            fprintf('No MIAWARA MA-data found\n');
            h2odata = xa;
        end
        clear xa
        
    case 'RS_THUN-MLS'
        xa = get_apriori_h2o_from_RS_THUN(mtime,[ptz(1,1) vmr(1)]); % RS-THUN + surface value (Zimm-meteo)
        xa_mls = get_h2o_apr_MLSclimat(mtime,Y.LATITUDE);
        
        if ~isempty(xa_mls)       
            h2odata = [xa(xa(:,1)>17000,:); xa_mls(xa_mls(:,1)<=10000,:)];
            
            clear xa_mls
        else
            fprintf('No MLS-data found\n');
            h2odata = xa;
        end
        clear xa
        
    case 'RS_THUN-USSTD'  
        xa = get_apriori_h2o_from_RS_THUN(mtime,[ptz(1,1) vmr(1)]); % RS-THUN + surface value (Zimm-meteo)
        h2odata  = load(fullfile('datafiles','H2O_USSTD.mat'));
        xa_usstd = h2odata;
        clear h2odata
        
        
        h2odata = [xa(xa(:,1)>17000,:); xa_usstd(xa_usstd(:,1)<=10000,:)];
        clear xa xa_usstd
        
    case 'MIA-MIA'        
        mia_t = get_miawara_level2_trop(mtime-12/24,mtime+12/24,3.1,1);
        if ~isempty(mia_t)
            [val,id]  = nanmin(abs(mia_t.time-mtime));
            fprintf('Closest tropo profile: %s\n',datestr(mia_t.time(id),31));
        else
            error('No MIAWARA-tropo data found!');           
        end

        mia_s = get_miawara_level2(datestr(floor(mtime),'yyyy-mm-dd'),datestr(ceil(mtime),'yyyy-mm-dd'),'224');

        if ~isempty(mia_s)
            fprintf('Closest strato profile: %s\n',datestr(mia_s.time,31));
        else
            mia_s = get_miawara_level2(datestr(floor(mtime),'yyyy-mm-dd'),datestr(ceil(mtime)+1,'yyyy-mm-dd'),'224');
            if ~isempty(mia_s)
                fprintf('Closest strato profile: %s\n',datestr(mia_s.time,31));
            else
                error('No MIAWARA-strato data found!');
            end
        end
        
        indt      = mia_t.p>17000;
        h2odata   = [[mia_t.p(indt) mia_t.vmr(indt,id)]; [mia_s.pressure mia_s.vmr]];
        R.iwv     = mia_t.iwv(id);
        clear indt mia
        
    case 'FTIR'
        t    = datenum([Y.YEAR Y.MONTH Y.DAY Y.HOUR Y.MINUTE Y.SECOND]);
        ftir = ftir_from_mysql(t-2/24,t+2/24);
        if isempty(ftir)
            Q = [];
            return
        end
        R.ftirtimes = ftir.time;
        
        snd  = get_apriori_h2o_from_RS_THUN(t,[ptz(1,1) vmr(1)]); % RS-THUN + surface value (Zimm-meteo)
        pf   = 10.^(interp1(ptz(:,3),log10(ptz(:,1)),ftir.alt));
        %pf   = 10.^(interp1(ptz(:,3),log10(ptz(:,1)),ftir.alt(ftir.alt<=10000)));
        h2odata = [snd(snd(:,1)>pf(1),:); [pf nanmean(ftir.h2o,1)'/1e6]; [snd(end,1) 1e-6]];
        %h2odata = [snd(snd(:,1)>pf(1),:); [pf nanmean(ftir.h2o(:,ftir.alt<=10000),1)'/1e6]; snd(snd(:,1)<pf(end),:)];
    otherwise
        error('R.H2O_apriori_profile was not understood!')
end

% add fixed-points
if isfield(R,'H2O_apc')&&~isempty(R.H2O_apc.fp) % add fp to a priori
    xa = h2odata;
    fp = R.H2O_apc.fp;
    
    % calculate p/z,T and vmr if necessary
    for i=1:size(fp,1)
        if isnan(fp(i,1))
            fp(i,1)    = 10.^(interp1(ptz(:,3),log10(ptz(:,1)),fp(i,3)));
        elseif isnan(fp(i,3))
            fp(i,3)   = interp1(log10(ptz(:,1)),ptz(:,3),log10(fp(i,1)));
        end
        if isnan(fp(i,2))
            fp(i,2)   = interp1(ptz(:,3),ptz(:,2),fp(i,3));
        else
            R.H2O_apc.T_ip = interp1(ptz(:,3),ptz(:,2),fp(i,3));
        end
        if isnan(fp(i,4))
            fp(i,4)  = rh2vmr(fp(i,1),fp(i,2),R.H2O_apc.cRH);
        end
        fprintf('Add fixed-point %.0fPa, %.2fK %.0fm, vmr: %.6f\n',fp(i,:));
    end

    R.H2O_apc.fp = fp;

    if ~isfield(R,'grid_corr')
        R.grid_corr = 1; % default
    end

    xf = fp(:,[1 4]);
    pf = fp(:,1);
    xf_min = nanmin(xf(:,1));
    xf_max = nanmax(xf(:,1));
    if xf_min>650e2
        xa     = [xa(1,:); xf; xa(xa(:,1)<xf_min,:)];
    else
        xa     = [xa(1,:); xa(xa(:,1)>xf_max&xa(:,1)<ptz_ground(1),:); xf; xa(xa(:,1)<xf_min,:)];
    end

    %%% Adapt p-grid!!!           
    p      = Q.P_GRID;
    if R.grid_corr>0
        for i=1:size(pf)
            p      = [p(p>pf(i)); pf(i); p(p<pf(i))];
            xid    = find(p==pf(i));

            if R.grid_corr==1
                % Check for too close grid-points
                lpm = interp1(ptz(:,1),ptz(:,3),p(xid-1:xid+1));

                if lpm(3)-lpm(2)<500
                    p = [p(1:xid); p(xid+2:end)];   
                end
                if lpm(2)-lpm(1)<500&&xid>2
                    p = [p(1:xid-2); p(xid:end)];
                end
            end
        end
    end
    Q.P_GRID = p;
    Q.ABS_SPECIES(1).GRIDS{1} = Q.P_GRID; %!!!
    h2odata = xa;
    clear xa
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
        SX_rel = xmlLoad (fullfile(curr_dir,'datafiles','v17_sx_apriori_H2O.xml'));
    case 'v18'
        SX_rel = xmlLoad(fullfile(curr_dir,'datafiles','v18_sx_apriori_H2O.xml'));
    case 'v22'
        SX_v22_abs = xmlLoad(fullfile(curr_dir,'datafiles','v22_sx_apriori_H2O_ppm.xml'));
        SX{1} = [3;0];
        SX{2}(:,1) = log10(Q.ABS_SPECIES(1).ATMDATA.GRID1);
        SX{2}(:,2) = interp1(SX_v22_abs(:,1),SX_v22_abs(:,2),Q.ABS_SPECIES(1).ATMDATA.GRID1);
        SX{2}(:,3) = 0.25;
        rel = SX{2}(:,2)./Q.ABS_SPECIES(1).ATMDATA.DATA;
        rel(rel<0.15 | rel>4) = 0.15;
        rel(rel>0.65) = 0.65;
        rel(isnan(rel)) = 0.65;
        rel(1) = rel(2);
        rel=smooth(rel,10);
        SX{2}(:,2) = rel.*Q.ABS_SPECIES(1).ATMDATA.DATA;
        clear rel SX_v22_abs
            
    case 'v22+tropo'
        if any(R.H2O_apc.pgrid==0)
            R.H2O_apc.pgrid(R.H2O_apc.pgrid==0) = ptz(1,1);
        end
        if length(R.H2O_apc.corrlength)==1
            R.H2O_apc.corrlength = ones(size(R.H2O_apc.pgrid))*R.H2O_apc.corrlength;
        end
        
        SX_v22_abs = xmlLoad(fullfile(curr_dir,'datafiles','v22_sx_apriori_H2O_ppm.xml'));
        
        
        SX_abs{2}(:,1) = Q.ABS_SPECIES(1).ATMDATA.GRID1;
        SX_abs{2}(:,2) = interp1(SX_v22_abs(:,1),SX_v22_abs(:,2),Q.ABS_SPECIES(1).ATMDATA.GRID1);
        SX_abs{2}(:,3) = 0.25;
        
        rel = SX_abs{2}(:,2)./Q.ABS_SPECIES(1).ATMDATA.DATA;
        rel(rel<0.15 | rel>4) = 0.15;
        rel(rel>0.65) = 0.65;
        rel(isnan(rel)) = 0.65;
        rel(1) = rel(2);
        rel=smooth(rel,10);
        SX_abs{2}(:,2) = rel;
        
        SX_rel{2}      = [[R.H2O_apc.pgrid(R.H2O_apc.pgrid>=500e2)' R.H2O_apc.std(R.H2O_apc.pgrid>=500e2)'...
            R.H2O_apc.corrlength(R.H2O_apc.pgrid>=500e2)'];[200e2 R.H2O_apc.std(end) R.H2O_apc.corrlength(end)]; SX_abs{2}(SX_abs{2}(:,1)<110e2,:)];
        SX_rel{2}(:,1) = log10(SX_rel{2}(:,1));
        SX_rel{1}      = [R.H2O_apc.corr_fun;R.H2O_apc.corr_cutoff];
        clear rel SX_v22_abs SX_abs
        
    case 'tropo'
        if any(R.H2O_apc.pgrid==0)
            R.H2O_apc.pgrid(R.H2O_apc.pgrid==0) = ptz(1,1);
        end
        fp = [];
        if R.H2O_apc.fix_surface
            fp = ptz(1,1);
        end    
        if ~isempty(R.H2O_apc.fp)
            fp = [fp; R.H2O_apc.fp(:,1)];
        end
        
        SX_rel = create_sa(R.H2O_apc.corr_fun,R.H2O_apc.corr_cutoff,R.H2O_apc.corrlength,R.H2O_apc.pgrid,R.H2O_apc.std,fp,Q.P_GRID);

    case 'file'
        SX_rel = xmlLoad('/tmp/SX.xml');
        
    otherwise
        error('R.H2O_apriori_covariance was not understood!')
end
if exist('SX_rel','var')
    switch Q.ABS_SPECIES(1).UNIT
        case 'vmr'
            SX{1}      = SX_rel{1};
            SX{2}(:,1) = log10(Q.ABS_SPECIES(1).ATMDATA.GRID1);
            SX{2}(:,2) = interp1(SX_rel{2}(:,1),SX_rel{2}(:,2),SX{2}(:,1)).*Q.ABS_SPECIES(1).ATMDATA.DATA;
            SX{2}(:,3) = interp1(SX_rel{2}(:,1),SX_rel{2}(:,3),SX{2}(:,1));
        case 'rel'
            SX         = SX_rel;
        case 'logrel'
            SX         = SX_rel;
    end
        
elseif ~exist('SX','var')
    error('No H2O covariance matrix defined!')
end
switch SX{1}(1)
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
Q.ABS_SPECIES(1).SX_rel = SX;
Q.ABS_SPECIES(1).SX = covmat1d_from_cfun( Q.ABS_SPECIES(1).GRIDS{1}, [power(10,SX{2}(:,1)) SX{2}(:,2)], corr_fun, [power(10,SX{2}(:,1)) SX{2}(:,3)], SX{1}(2), @log10);
clear SX_rel SX h2odata


clear ptz

 %- Other species
% Oxygen
O2 = xmlLoad(fullfile(curr_dir,'datafiles','apriori_O2.xml'));
Q.ABS_SPECIES(2).TAG       = {'O2-PWR93'};
Q.ABS_SPECIES(2).RETRIEVE  = false;
Q.ABS_SPECIES(2).ATMDATA.TYPE       = 'atmdata';
Q.ABS_SPECIES(2).ATMDATA.NAME       = 'O2';
Q.ABS_SPECIES(2).ATMDATA.SOURCE     = 'unknown';
Q.ABS_SPECIES(2).ATMDATA.DIM        = 1;
Q.ABS_SPECIES(2).ATMDATA.DATA       = O2(:,2);
Q.ABS_SPECIES(2).ATMDATA.DATA_NAME  = 'Volume mixing ratio';
Q.ABS_SPECIES(2).ATMDATA.DATA_UNIT  = '-';
Q.ABS_SPECIES(2).ATMDATA.GRID1      = O2(:,1);
Q.ABS_SPECIES(2).ATMDATA.GRID1_NAME = 'Pressure';
Q.ABS_SPECIES(2).ATMDATA.GRID1_UNIT = 'Pa';

% Nitrogen
N2 = xmlLoad(fullfile(curr_dir,'datafiles','apriori_N2.xml'));
Q.ABS_SPECIES(3).TAG       = {'N2-SelfContStandardType'};
Q.ABS_SPECIES(3).RETRIEVE  = false;
Q.ABS_SPECIES(3).ATMDATA.TYPE       = 'atmdata';
Q.ABS_SPECIES(3).ATMDATA.NAME       = 'N2';
Q.ABS_SPECIES(3).ATMDATA.SOURCE     = 'unknown';
Q.ABS_SPECIES(3).ATMDATA.DIM        = 1;
Q.ABS_SPECIES(3).ATMDATA.DATA       = N2(:,2);
Q.ABS_SPECIES(3).ATMDATA.DATA_NAME  = 'Volume mixing ratio';
Q.ABS_SPECIES(3).ATMDATA.DATA_UNIT  = '-';
Q.ABS_SPECIES(3).ATMDATA.GRID1      = N2(:,1);
Q.ABS_SPECIES(3).ATMDATA.GRID1_NAME = 'Pressure';
Q.ABS_SPECIES(3).ATMDATA.GRID1_UNIT = 'Pa';

% CO2
CO2 = xmlLoad(fullfile(curr_dir,'datafiles','apriori_CO2.xml'));
Q.ABS_SPECIES(4).TAG       = {'CO2-SelfContPWR93'};
Q.ABS_SPECIES(4).RETRIEVE  = false;
Q.ABS_SPECIES(4).ATMDATA.TYPE       = 'atmdata';
Q.ABS_SPECIES(4).ATMDATA.NAME       = 'CO2';
Q.ABS_SPECIES(4).ATMDATA.SOURCE     = 'unknown';
Q.ABS_SPECIES(4).ATMDATA.DIM        = 1;
Q.ABS_SPECIES(4).ATMDATA.DATA       = CO2(:,2);
Q.ABS_SPECIES(4).ATMDATA.DATA_NAME  = 'Volume mixing ratio';
Q.ABS_SPECIES(4).ATMDATA.DATA_UNIT  = '-';
Q.ABS_SPECIES(4).ATMDATA.GRID1      = CO2(:,1);
Q.ABS_SPECIES(4).ATMDATA.GRID1_NAME = 'Pressure';
Q.ABS_SPECIES(4).ATMDATA.GRID1_UNIT = 'Pa';

%- Define L2 structure (beside retrieval quantities below)
Q.L2_EXTRA = {'dx','cost','e','eo','es','yf','A','S','So','Ss','mresp','date','xa','y','bl','tnoise','ptz','G','J'};
