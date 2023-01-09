% function Y = msm2Y(msm)
%
% This function creates the Y-Structure, required by
% Qpack2, based on information of the measurements (msm)
%
% DS, 2012-01-09
%
function Y = msm2Y(msm)
    
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
Y.Z_PLATFORM = 16000;            % Platform altitude [m]
Y.ZA         = 0;                % Zenith angle of the measurement [deg]
Y.HSE_P      = 10000;            % Hydrostatic equilibrium reference pressure [Pa]
Y.HSE_Z      = 16000;            % Hydrostatic equilibrium reference height [m]
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
