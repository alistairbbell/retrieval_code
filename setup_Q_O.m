%%%%%%%%% define forward model and retrieval settings %%%%%%%%%
%[Q, O] = setup_Q_O
function [Q, O] = setup_Q_O(Y,msm, extra_info)

%=== Directories
top_dir = extra_info.top_dir;
arts_includes  = atmlab( 'ARTS_INCLUDES' );
atmlab( 'VERBOSITY', 3);
atmlab( 'FMODEL_VERBOSITY', 0);

R.PTZ_profile = 'ECMWF';
ptz = get_ptz_ECMWF_tub_temporal_interpolation(msm, extra_info,Y);
disp('ptz(1,:)')
disp(ptz(1,:))

p_smooth = 10.^ linspace(5.0,-2.25,145);
t_smooth = interp1(ptz(:,1),smooth(ptz(:,2),4),  p_smooth);
z_smooth = interp1(ptz(:,1), ptz(:,3),  p_smooth);
ptz_smooth = [p_smooth; t_smooth; z_smooth]';

%%%%%%%%% initialize Q
Q = qarts;
Q = createQ(Y,R,ptz);

%- General
Q.INCLUDES            = { fullfile( arts_includes, 'general.arts'),fullfile( arts_includes, 'continua.arts' ) ,...
fullfile( arts_includes, 'agendas.arts' ), fullfile( 'ARTS_INCLUDES', 'planet_earth.arts' ) };
Q.ATMOSPHERE_DIM      = 1;
Q.STOKES_DIM          = 1;
Q.J_DO                = true;
Q.CLOUDBOX_DO         = false;
Q.INPUT_FILE_FORMAT   = 'double';
Q.OUTPUT_FILE_FORMAT  = 'ascii';

%- Radiative transfer
%
Q.IY_UNIT              = 'RJBT'; 
Q.YCALC_WSMS          = { 'yCalc' };
%
Q.PPATH_LMAX          = 250; % max distance between points
%Q.PPATH_STEP_AGENDA   = { 'ppath_stepGeometric' };  % geom. Propagation, 
                                                    % refraction neglected
Q.PPATH_AGENDA               = { 'ppath_agenda__FollowSensorLosPath'   };
Q.PPATH_STEP_AGENDA          = { 'ppath_step_agenda__GeometricPath'    };
Q.IY_SPACE_AGENDA            = { 'iy_space_agenda__CosmicBackground'   };
Q.IY_SURFACE_AGENDA          = { 'iy_surface_agenda__UseSurfaceRtprop' };
Q.IY_MAIN_AGENDA             = { 'iy_main_agenda__Emission'            };

%- Surface
Q.Z_SURFACE          = 15e3;      % Just a dummy value. A 15 km observation 
                                  % altitude is assumed here

%- Absorption
%
Q.ABS_LINES           = [top_dir,'/nedoluha1995_hyperfine.xml'];% fullfile( atmlab_example_data, 'o3line111ghz' )
Q.ABS_LINES_FORMAT    = 'Arts';
Q.ABS_LINES_FORMAT    = 'ARTSCAT';
Q.ABSORPTION          = 'OnTheFly';
Q.ABS_NLS             = [];

z_toa                 = 150e3;
Q.SENSOR_LOS          = Y.ZA;
Q.SENSOR_POS          = Y.Z_PLATFORM;
Q.SENSOR_DO           = true;
% Frequency, spectrometer and pencil beam antenna

%Define frequency grid for the forward model (more dense in centre, less
%dene at wings
fgrid = expStepApprox(extra_info.f_cen,(extra_info.f_cen+1e6+extra_info.bw_fm/2),5e4,1200 );
fgrid_firstDiff = fgrid(2)-fgrid(1);
fgrid_2 = expStepApprox(extra_info.f_cen-fgrid_firstDiff,(extra_info.f_cen-extra_info.bw_fm/2), 5e4,1200);
Q.F_GRID              = [fgrid_2; fgrid];
%Define Backend properties (backend do must be set to true)

% The hypothetical spectrometer has rectangular response functions
ch_resp = read_datafile(fullfile(top_dir,'channelresponse.aa'),'matrix');
H                     = qartsSensor;  % initialize sensor
H.SENSOR_NORM         = true;
H.F_BACKEND   = Y.F;%[ min(Q.F_GRID)+df : df : max(Q.F_GRID)-df )]';
H.BACKEND_DO = true;

B.name                = 'Spectrometer channel response function';
B.gridnames           = {'Frequency'} ;
B.grids               =  {ch_resp(:,1)'};
B.dataname            = 'Response';
B.data                = ch_resp(:,2)';
%
H.BACKEND_CHANNEL_RESPONSE{1} = B;
clear B ch_resp

% Acqiris FFT channel response
% ch_resp = read_datafile(fullfile(curr_dir,'datafiles','FFT_ACQ1_channelresponse.aa'),'matrix');
% Q.SENSOR_RESPONSE.BACKEND_DO = 1;
% Q.SENSOR_RESPONSE.MIXER_DO = 0;
% Q.SENSOR_RESPONSE.F_BACKEND  = Y.F;
% Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.name      = 'Acqiris FFT channel response';
% Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.gridnames = {'Frequency'};
% Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.grids     = {ch_resp(:,1)'};
% Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.dataname  = 'Response';
% Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.data      = ch_resp(:,2)';
% clear ch_resp

%
Q.SENSOR_DO           = true; %to include sensor settings
Q.SENSOR_RESPONSE     = H;
Q.ANTENNA_DIM         = 1;
Q.MBLOCK_DLOS_GRID    = 0;

%- Correlation of thermal noise (not correlated)
Q.TNOISE_C = covmat1d_from_cfun(Q.SENSOR_RESPONSE.F_BACKEND,[],'drc',0,0);  %no correlation, diagonal in Y.TNOISE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                Pressure & Height Variables                         %%%
%- Pressure grid (not the retrieval grid!)
Q.P_GRID   = power(10,linspace(4.25,-2,100))';
Q.HSE.ON       = false;

if isnan(mean(ptz(:,2), 'omitnan'))
    disp('No PTZ data found!')
    return
end

Q.VMR_NEGATIVE_OK = false; %avoid negative values of water vapour in forward model
Q.T.RETRIEVE           = false;
Q.T.ATMDATA.TYPE       = 'atmdata';
Q.T.ATMDATA.NAME       = 'Temperature';
Q.T.ATMDATA.SOURCE     = 'MLS climatology';
Q.T.ATMDATA.DIM        = 1;
Q.T.ATMDATA.DATA       = ptz_smooth(:,2);
Q.T.ATMDATA.DATA_NAME  = 'Temperature';
Q.T.ATMDATA.DATA_UNIT  = 'K';
Q.T.ATMDATA.GRID1      = ptz_smooth(:,1);
Q.T.ATMDATA.GRID1_NAME = 'Pressure';
Q.T.ATMDATA.GRID1_UNIT = 'Pa';
Q.Z.ATMDATA.TYPE       = 'atmdata';
Q.Z.ATMDATA.NAME       = 'Altitude';
Q.Z.ATMDATA.SOURCE     = 'MLS';
Q.Z.ATMDATA.DIM        = 1;
Q.Z.ATMDATA.DATA       = ptz_smooth(:,3);
Q.Z.ATMDATA.DATA_NAME  = 'Altitude';
Q.Z.ATMDATA.DATA_UNIT  = 'm';
Q.Z.ATMDATA.GRID1      = ptz_smooth(:,1);
Q.Z.ATMDATA.GRID1_NAME = 'Pressure';
Q.Z.ATMDATA.GRID1_UNIT = 'Pa';
clear ptz
clear ptz_smooth

%%% a priori profile
%Q.ABS_SPECIES(1).ATMDATA  = gf_load([top_dir,'/h2o_climatology_mls22.mat']); does not yet work
%for i = 1:length(Y.YEAR)

disp('DATETIME A PRIORI')
disp(datetime(Y.YEAR(1), Y.MONTH(1), Y.DAY(1), Y.HOUR(1), Y.MINUTE(1), Y.SECOND(1)))
disp('getting h20')

%get h2o profile from zimmerwald (Bern) climatology
h2o=get_apriori_h2o_from_ecmwf_zimmerwald(datetime(Y.YEAR(1), Y.MONTH(1), Y.DAY(1), Y.HOUR(1), Y.MINUTE(1), Y.SECOND(1)),extra_info.a_priori_filename );

[pn,ind]=sort(-h2o(:,1));
h2on=h2o-h2o;
h2on(:,1)=h2o(ind,1);
h2on(:,2)=h2o(ind,2);
h2o=h2on; 

%- Species
% Water vapor
Q.ABS_SPECIES(1).TAG      = {'H2O'};
Q.ABS_SPECIES(1).RETRIEVE = true;
Q.ABS_SPECIES(1).L2       = true;
Q.ABS_SPECIES(1).GRIDS{1} = power(10,linspace(4,-1,80))'; %power(10,linspace(4,0,50))';somewhat denser than old grid
Q.ABS_SPECIES(1).GRIDS{2} = [Y.LATITUDE]; %lat grid
Q.ABS_SPECIES(1).GRIDS{3} = []; %lon grid, for zonal mean, lon=0
Q.ABS_SPECIES(1).MINMAX   = 1e-8;
Q.ABS_SPECIES(1).UNIT     = 'vmr';

Q.ABS_SPECIES(1).ATMDATA.TYPE       = 'atmdata';
Q.ABS_SPECIES(1).ATMDATA.NAME       = {'H2O'};
Q.ABS_SPECIES(1).ATMDATA.SOURCE     = 'ECMWF climatology';
Q.ABS_SPECIES(1).ATMDATA.DIM        = 1;
Q.ABS_SPECIES(1).ATMDATA.DATA       = h2o(:,2);
Q.ABS_SPECIES(1).ATMDATA.DATA_NAME  = '';
Q.ABS_SPECIES(1).ATMDATA.DATA_UNIT  = 'ppm';      
Q.ABS_SPECIES(1).ATMDATA.GRID1   = h2o(:,1);
Q.ABS_SPECIES(1).ATMDATA.GRID1_NAME = 'Pressure';
Q.ABS_SPECIES(1).ATMDATA.GRID1_UNIT = 'Pa';

% apriori profile now with descending pressure 1000, ...., 0.01!!!!!!



% read apriori covariance matrix relative and absolute
sx_ppm=read_datafile([top_dir,'/sx.apriori.H2O.ppm.aa'],'matrix');
sx_ppm(:,2) = ones(1,85)*4e-5;

% convert sigma from ppm to fraction of a priori profile, on a priori grid
sx_rel(:,1)=Q.ABS_SPECIES(1).ATMDATA.GRID1;
sx_rel(:,2)=interp1(log10(sx_ppm(:,1)),sx_ppm(:,2),log10(Q.ABS_SPECIES(1).ATMDATA.GRID1))./h2o(:,2);

% constrain the relative sigma to 30-60%
sx_rel(sx_rel(:,2)>0.6 ,2)=0.6;
sx_rel(sx_rel(:,2)<0.3 ,2)=0.3;

sx_perc = ones(length(sx_rel(:,1)),1);
sx_perc(1:618,:) = 0.2;
sx_perc(618:943,:) = linspace(0.2,0.4,(943-617));
sx_perc(943:end) = 0.4;

% % smooth the relative sigma profile
sx_rel(:,2)= sx_perc;
sx_rel(:,2)=smooth(sx_rel(:,2),10);

% %%sx_rel=read_datafile([top_dir,'/sx.apriori.H2O.rel.aa'],'matrix');    sx
% %%with fixed rel
% 
sx_abs=[sx_rel(:,1),sx_rel(:,2).*Q.ABS_SPECIES(1).ATMDATA.DATA];        % must match Q.ABS_SPECIES.UNIT

plot_a_prioiri_errs(sx_abs, extra_info.plot_dir)
%sx_abs(:,1)=Q.ABS_SPECIES(1).ATMDATA.GRID1;
%sx_abs(:,2)=interp1(log10(sx_ppm(:,1)),sx_ppm(:,2),log10(Q.ABS_SPECIES(1).ATMDATA.GRID1));
corr_length=0.25;
%corr_len_2d = [sx_abs(:,1)  linspace(0.01,corr_length, length(sx_abs(:,1)))'];
%corr_len_2d = [sx_abs(:,1)  ones(length(sx_abs(:,1)),1)*corr_length];
%corr_len_2d(1:943,2) =  linspace(0.001,corr_length, 943)';
corr_len_2d(1:943,2) =  0.1;
plot_corr_len(corr_len_2d, extra_info.plot_dir)

%disp(Q.ABS_SPECIES(1).ATMDATA.GRID1)
disp('size(Q.ABS_SPECIES(1).GRIDS{1})')
disp(size(Q.ABS_SPECIES(1).GRIDS{1}))

% write it to file
Q.ABS_SPECIES(1).SX = covmat1d_from_cfun( Q.ABS_SPECIES(1).GRIDS{1},sx_abs, 'exp',corr_length, 0, @log10);
% Polyfit
% A polynomial of order 2 is used for "baseline fit".

Q.POLYFIT.RETRIEVE        = true;
Q.POLYFIT.ORDER           = 4;
Q.POLYFIT.L2              = true;
Q.POLYFIT.SX0             = 1^2;
Q.POLYFIT.SX1             = 0.1^2;
Q.POLYFIT.SX2             = 0.1^2;
Q.POLYFIT.SX3             = 0.1^2;
Q.POLYFIT.SX4             = 0.1^2;

%- Define L2 structure (beside retrieval quantities below)
Q.L2_EXTRA = {'dx','cost','e','eo','es','A','S','yf','So','Ss','G','J','mresp','date','xa','y','bl','tnoise','ptz'};
%
%- Frequency (shift retrieval here, shift+stretch also possible)
%
Q.FSHIFTFIT.RETRIEVE = true;
Q.FSHIFTFIT.DF = 5e3;
Q.FSHIFTFIT.SX = 200e2^2;
Q.FSHIFTFIT.L2 = true;
Q.SINEFIT.RETRIEVE = false;
Q.SINEFIT.PERIODS = [ 80e6 ]';%[ 5e6 80e6 ]
Q.SINEFIT.L2 = true;
Q.SINEFIT.SX1 = 5e6;% a priori error covariance matrix for coefficients of first period length
%Q.SINEFIT.SX2 = 10e6;

%- Frequency stretch
Q.FSTRETCHFIT.RETRIEVE = false;% Set to true to activate retrieval

disp('size(Q)')
disp(size(Q))

%%% Define O
O = qp2_l2( Q );  % This ensures that OEM returns the varibles needed
O.linear = false;
O.stop_dx = 1e-3;
O.yf = true;
%O.jexact = true;


function output = expStepApprox(A,B,f1,N)
output = zeros(N,1);
step_in = log2(f1);
if B>A
    step_fin = log2(B-A);
    steps = linspace(step_in,step_fin,N);
elseif A>B
    step_fin = log2(A-B);
    steps = linspace(step_fin, step_in,N);
else
    disp('A and B are the same! No frequency steps created!')
end
for i = 1:N
    if B>A
        output(i) = (A -f1 + 2^(steps(i)));
    else
        output(i) = (A + f1 - 2^(steps(i)));
    end
end    

end



end
