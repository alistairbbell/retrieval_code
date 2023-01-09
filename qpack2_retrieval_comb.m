% function L2 = qpack2_retrieval(M,R,QQ)
%
% This function does a Qpack2 retrieval, based on information
% about the measurement (M) and on the information about the
% retrieval (R). The QQ-Structure contains additional information
% for the Q-Structure. The information in the QQ-Structure will be
% replaced into the Q-Structure after the Q-Structure was created.
%
function L2 = qpack2_retrieval_comb(M,R,QQ)
%%% Init atmlab
run( '/opt/arts/atmlab/atmlab/atmlab_init' );

L2 = [];

%%% Load the measurements
disp(' ');
disp('Loading the measurements...')
%msm = loadmsm(M);
load(M.filename);
if isempty(msm)||isempty(msm.y)
    disp('No measurements!!')
    return
end

%%% Create the Y-Structure required by Qpack2
disp('Creating the Y-structure...')
Y = msm2Y_tr(msm);

%%% Create the Q-Structure required by Qpack2
disp('Creating the Q-structure...')
[Q,R] = createQ_tr(Y,R);

if isempty(Q)
    disp('Could not create Q structure!')
    return
end

if exist('QQ','var') && ~isempty(QQ)
    QQnames = fieldnames(QQ);
    for i = 1:size(QQnames,1)
        eval(['Q.' QQnames{i} ' = QQ.' QQnames{i} ';']);
    end
end
Q.SENSOR_RESPONSE.MIXER_DO   = 0;
Q.SENSOR_RESPONSE.BACKEND_DO = 0;

%%% Check that all frequencies are OK
Y = rmfield(Y,'MIN_TIME');
Y = rmfield(Y,'MAX_TIME');
if Q.SENSOR_RESPONSE.MIXER_DO
    Y.F = Q.SENSOR_RESPONSE.F_BACKEND;
end
%disp('Checking for consistency...')
%if ~qp2_check_f( Q, Y, 1e3 );
    %disp( 'Some mismatch between Q.F_BACKEND and frequencies of spectra.' );
    %return
%end


%%% Define OEM variables
O = qp2_l2( Q );  % This ensures that OEM returns the variables needed
if ~R.linear
    O.linear           = false;
    O.ga_start         = 100;
    O.ga_factor_ok     = 10;
    O.ga_factor_not_ok = 10;
    O.stop_dx          = 0.1;
    O.ga_max           = 1e5;
    O.maxiter          = 20;
    O.Xiter            = true;
else
    O.linear           = true;
end        

%%% Make (1st) inversion
disp('Doing inversion...')
L2 = qpack2( Q, O, Y );
if isempty(L2)
    disp('Could not do the inversion!')
    return
end

%%% Reset measurement noise
% For total power spectra, the measurement noise (sigma) is determined 
% using the difference between the forward model of a first guess retrieval
% and the measured spectrum 


%%% Process L2 (+characterisation)
n    = length(L2.species1_x);

% add setups
L2.M = M;
L2.R = R;
L2.O = O;
L2.Q = Q;
L2.Y = Y;
%L2.Q.SENSOR_RESPONSE.F_BACKEND = L2.Q.SENSOR_RESPONSE.F_BACKEND + L2.Q.SENSOR_RESPONSE.LO;

% Se, Sx
[L2.char.Se,L2.char.Sx]  = get_Se_Sx_from_L2(L2,1);

%%% stuff for rel/logrel only
if strcmp(Q.ABS_SPECIES(1).UNIT,'rel')||strcmp(Q.ABS_SPECIES(1).UNIT,'logrel') 
    
    % Change L2 from rel/logrel
    L2.char.mrunsc = L2.species1_mr;
    L2.char.Aunsc  = L2.species1_A;
    L2 = qp2_rel2vmr(L2);  
    
    % Recalculate J,G,A with scaled J
    [L2.char.Js,L2.char.Gs,As] = scale_ret(full(L2.J),L2.species1_x./L2.species1_xa,L2.char.Se,L2.char.Sx,'scale');
    L2.species1_Arel  = As(1:n,1:n);
    L2.species1_mrrel = sum( L2.species1_Arel' );

    % relative errors
    Sobs          = L2.char.Gs*L2.char.Se*L2.char.Gs';
    L2.char.eorel = sqrt(diag(Sobs(1:n,1:n)));
    Ai            = As-eye(size(As));
    Ssmo          = Ai*L2.char.Sx*Ai';
    L2.char.esrel = sqrt(diag(Ssmo(1:n,1:n)));
    L2.char.erel  = sqrt(diag(Sobs(1:n,1:n)+Ssmo(1:n,1:n)));
    
    A             = L2.species1_Arel;
    mr            = L2.species1_mrrel;    
else
    A             = L2.species1_A;
    mr            = L2.species1_mr;
end

% z-grid
L2.species1_z = interp1(Q.Z_ATMDATA.GRID1,Q.Z_ATMDATA.DATA,L2.species1_p);

% level 1 data
L2.level1          = msm;
L2.level1.sigma    = Y.TNOISE;
L2.min_time        = msm.min_time;
L2.max_time        = msm.max_time;
L2.min_time_supposed = msm.startdate;
L2.max_time_supposed = msm.enddate;

disp('Doing characterisation')
% information content and characterisation
if length(L2.f)<1000 % for larger Se it takes too long...
    Ks             = sqrtm(L2.char.Se)\L2.char.Js(:,1:n)*sqrtm(L2.char.Sx(1:n,1:n));
    SKs            = svd(Ks);
    L2.char.rank   = length(find(SKs>1));     % rank
end
L2.char.ds     = trace(A);    % degrees of freedom
L2.char.H      = -1/2*log(prod(eig(eye(n)-A))); % Shannon information content
L2.char.ap     = get_avk_params(A,mr,L2.species1_p,L2.species1_z);
L2.char.offset = L2.polyfit0_x;
L2.char.slope  = L2.polyfit1_x;

% meteo data
[L2.R.cloud,L2.R.prcp]=get_meteo_tags(msm.min_time,msm.max_time);
end

function Y = msm2Y_tr(msm)
    
%%% Create the Y structure requiered by Qpack2
Y            = qp2_y;
Y.YEAR       = year(msm.time);
Y.MONTH      = month(msm.time);
Y.DAY        = day(msm.time);
Y.HOUR       = hour(msm.time);
Y.MINUTE     = minute(msm.time);
Y.SECOND     = second(msm.time);
Y.LATITUDE   = msm.latitude;     % Instrument latitude [deg]
Y.LONGITUDE  = msm.longitude;    % Instrument longitude [deg]
Y.Z_PLATFORM = msm.altitude;     % Platform altitude [m]
Y.ZA         = msm.za;           % Zenith angle of the measurement [deg]
Y.HSE_P      = msm.hse_p;        % Hydrostatic equilibrium reference pressure [Pa]
Y.HSE_Z      = msm.hse_z;        % Hydrostatic equilibrium reference height [m]
Y.F          = msm.f;            % Frequency vector [Hz]
Y.Y          = msm.y;            % Measurement vector, Brightness temperature [K]
Y.TNOISE     = msm.sigma;        % Measurement noise [K], i.e std(Y.Y)

%%% The following two variables are not "official" variables for QPack2 and
%%% must be removed before the retrieval. They are only used for internal
%%% purposes.
Y.MIN_TIME   = msm.min_time; 
Y.MAX_TIME   = msm.max_time;

%%% Adjust the measurement and frequency vector orientation, if necessary.
if size(Y.Y,1)<size(Y.Y,2); Y.Y = Y.Y'; end
if size(Y.F,1)<size(Y.F,2); Y.F = Y.F'; end
end
