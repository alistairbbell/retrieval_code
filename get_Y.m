%%%%%%%%% Import measurement data and set up Y
% FORMAT   [Y, msm] = get_Y(start_num,narg2,bw)
%        
% OUT   Y        
%       msm      level1 data
% IN    start_num
% 	narg2	 either sigma or end_num
%	bw	 bandwidth of measurement

function [Y, msm] = get_Y(start_num,narg2,bw)

disp('-> loading measurement')
bandwidth=bw;
msm.y=[];
Y            = qp2_y;               % initialize

mission = get_mission_from_mysql(start_num);
if ~isempty(mission)
    
if narg2<1
	sigma = narg2;
	
	if start_num>datenum(2010,12,15)
		disp('-> OMT measurement')
		msm = integrate_level1_dual_combo_sigma(datestr(start_num,31),sigma,mission);
	else
		msm=integrate_level1_acq_from_mysql(datestr(start_num,31),sigma,mission);
	end
else
    end_num = narg2;
    if start_num>datenum(2010,12,15)
		disp('-> OMT measurement')
		msm = integrate_level1_dual_combo(datestr(start_num,31),datestr(end_num,31),mission);
	else
		msm=get_level1_acq_from_mysql_certain_time_periode(datestr(start_num,31),datestr(end_num,31),mission);
    end
end
end

if ~isempty(msm.y)==1
    frequency_ind = (find(msm.f>=22.23508e9-bandwidth/2 & msm.f<=22.23508e9+bandwidth/2));
    msm.y = msm.y(frequency_ind);
    msm.f = msm.f(frequency_ind);
else
   disp('no level1 data found, inversion aborted!')
   msm=[];
   return
end


%%% set up Y as required by qpack2:

    Y.YEAR       = year(msm.time);
    Y.MONTH      = month(msm.time);
    Y.DAY        = day(msm.time);
    Y.HOUR       = hour(msm.time);
    Y.MINUTE     = minute(msm.time);
    Y.SECOND     = second(msm.time);
    Y.LATITUDE   = msm.latitude;
    Y.LONGITUDE  = msm.longitude;
    Y.HSE_P      = 100*nanmean(mysql(sprintf('select pressure from meteo where time between "%s" and "%s";', datestr(msm.time-0.1,31),datestr(msm.time+0.1,31)),'mat'));
    delta=1;
    while isnan(Y.HSE_P) %if pressure is missing, increase time interval gradually
        Y.HSE_P      = 100*nanmean(mysql(sprintf('select pressure from meteo where time between "%s" and "%s";', datestr(msm.time-delta,31),datestr(msm.time+delta,31)),'mat'));
        delta=delta+1;
    end
    Y.HSE_Z      = msm.altitude;
    
    Y.Z_PLATFORM = 16000; % Platform altitude
    Y.ZA         = 0;     % Zenith angle of the measurement
    Y.F = msm.f-100e3;
    Y.Y = msm.y';
    
    %%% limit noise to noise >=0.01K!
    if msm.sigma>0.01
        Y.TNOISE = msm.sigma; % Measurement noise, i.e std(Y.Y)
    else
        Y.TNOISE = 0.01;
    end
