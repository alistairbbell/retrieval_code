%--------------------------------------------------------------------------
% Get radiosounding data from mysql
%  additionally, vmr(volume mixing ratio in [ ]), rho (density in [kgm-3])
%  and IWV (column density in [mm]) are calculated
%
% rsd = get_snd_from_mysql(startdate,enddate,stnid,prec,corr)
%
%  Input:   startdate,endate            timeframe (in MATLAB time format, datenum(:))
%           stnid                       stationid (Payerne=102)
%           prec                        flag to set timeframe behaviour:
%                                       0: make sure to include only entire profiles of nighttime soundings
%                                       1: take the timeframe as defined with startdate,enddate (and risk to only get partial 
%                                          profiles of nighttime soundings)
%           corr                        0: do not apply Miloshevich RS correction
%                                       1: apply day and nightime correction
%                                       2: apply only daytime correction
%                                       3: apply only daytime correction and not below clouds
%
%  Output:  rsd                         Structure containing vectors and cells with sounding profile data, 
%                                        one cell per sounding
%                                       starttime     vector of launching times (in MATLAB-format)
%                                       pres{i}       pressure  [Pa] of sounding i
%                                       alt{i}        altitude  [m]
%                                       T{i}          temperature [K]
%                                       RH{i}         relative Humidity [%]
%                                       rho{i}        water vapour density [kgm-3]
%                                       vmr{i}        water vapour volume mixing ratio [ ]
%                                       rho_air{i}    air density [kgm-3]
%                                       iwv           vector of water vapour column densities [mm(=kgm-2)]
%
%                                       .._raw{i}     uncorrected profile data (only if RH-correction was applied, 
%                                                     d.h. if corr~0)
%                                       iwv_raw       iwv of uncorrected profile [mm(=kgm-2)]
%                                       wdir          wind direction [°]
%                                       wv            wind velocity [ms-1]
%                                       
% (c) Rene Bleisch, rene.bleisch@iap.unibe.ch 2011
%--------------------------------------------------------------------------
function rsd = get_snd_from_mysql(startdate,enddate,stnid,prec,corr)
if nargin<5
    corr = 0;
end

% Some useful constants
Ma      = .028966;   % Molar mass of dry air        [kg mol-1]
R       = 8.314472;  % Universal gas constant       [J mol-1 K-1]
Rw      = R/Ma;      % Gas constant of dry air      [J kg-1 K-1]

% Known bad Payerne profiles (unrealistic RH-values)
badsnds = [datenum(2011,1,16,11,01,00),...
    datenum(2011,2,12,11,03,00)];

%=== Get raw data from mysql ==============================================
if prec == 0 
%   As there unfortunately exists no profile_ID, soundings only can be hold apart by using time information
%   prec=0 can be used to avoid getting only partial profiles:
%   Example:   get_snd_from_mysql(datenum(2009,1,1),datenum(2009,1,2),...) would normally deliver a part of the profile from the 
%   31.12.2008 ~23:00 sounding (only values from 1.1.2009 0:00, the lowermost part is missing) and a part of the profile 
%   from the 2.1.2009 ~23:00 sounding(only the lowermost part).
%   Using prec=0 makes sure to include the entire profiles of these soundings. 

    startdate =floor(startdate)-1+22/24; % usually soundings are started ~23:00...
    enddate   =floor(enddate)+1+2/24;    % ... and finished at least ~1:00
end

% change data format
startdate = datestr(startdate,31);
enddate   = datestr(enddate,31);

mysql('use RADIOSONDE;');

% SELECT string
ystr        = 'year(datetime) as year';
mstr        = 'month(datetime) as month';
dstr        = 'dayofmonth(datetime) as day';
tstr        = 'hour(datetime) as hour, minute(datetime) as minute,second(datetime) as second';
selstr      = sprintf('%s, %s, %s, %s, height, pressure, temperature, rel_humidity, wind_dirn, wind_speed',ystr,mstr,dstr,tstr);


% WHERE string 
 % stationid and datetime
wstr        = sprintf('stationid = "%d" AND datetime between "%s" AND "%s"',stnid,startdate,enddate);

% assemble and execute mysql-command 
cmd         = sprintf('SELECT %s FROM radiosonde_measurement WHERE %s ORDER by datetime',selstr,wstr);
data        = mysql(cmd,'mat');

% get time and height
time        = nan(length(data(:,1)),1);
for i = 1:length(data(:,1))
    time(i)        = datenum(data(i,1:6));
end
height      = data(:,7);
ndata       = length(time);

if ndata == 0
    rsd     = [];
    return
end

%=== Split raw data into separate profiles ====================================
% (Necessary, because there exists no key in mysql-db which clearly separates each profile :-( )

nprof       = 1; % # of profiles
sind        = 1; % indices, at which each sounding starts
eind        = [];% indices, at which each sounding ends

% assumes a new sounding, if time(i)-time(i-1) is greater than 1h
for i = 2:ndata
    if time(i)-time(i-1)>0.5/24
        nprof = nprof + 1;
        sind = [sind; i];
        eind = [eind; i-1];
    end
end
eind = [eind; ndata];

times = nan(nprof,1);
iwv   = nan(nprof,1);

% make one cell per sounding (a cell variable can consist of cells of different size 
%  (each profile has a different length)
q = 1;
for p = 1:nprof
    % Check for bad snds
    if stnid==102
        if any(badsnds==time(sind(p)))
            if nprof>1
                continue
            else
                rsd = [];
                return
            end
        end
    end
    times(q)= time(sind(p));
    alt{q} = height(sind(p):eind(p));
    pres{q}= data(sind(p):eind(p),8)*100;
    T{q}   = data(sind(p):eind(p),9);
    RH{q}  = data(sind(p):eind(p),10);
    
    % find RH<=1 and set to NaN (untrusty)
    RH{q}(RH{q}<=1) = NaN;
    
    vmr{q}  = humidity_transform(RH{q},'RH','VMR',T{q},pres{q});
    rho{q}  = humidity_transform(RH{q},'RH','HUM',T{q},pres{q});
    rho_air{q} = pres{q}./(Rw*T{q});

    % calculate iwv
    ind     = find(~isnan(rho{q})); % humidity data is available only up to 10-12 km
    if length(ind)>1
        iwv(q)  = trapz(alt{q}(ind),rho{q}(ind)); 
    end
    
    % wind
    wdir{p} = data(sind(p):eind(p),11);
    wv{p}   = data(sind(p):eind(p),12);
    
    q       = q+1;
        
end

if q-1~=nprof    
    nprof     = q-1;
    times     = times(1:nprof);
    iwv       = iwv(1:nprof);  
end

rsd.alt       = alt;
rsd.starttime = times;
rsd.p         = pres;
rsd.T         = T;

%=== RS92_correction according to Miloshevich 2009 ==============================
if corr ~=0
    % Get lat&lon
    cmd   = sprintf('SELECT latitude,longitude from startwaveh2o.station WHERE stationid="%d"',stnid);
    coord = mysql(cmd,'mat');

    iwv_raw = iwv;

    for p = 1:nprof
        RH_raw{p} = RH{p};
        RH{p}     = correct_rs92(times(p),coord(1),coord(2),pres{p}/100,RH{p},corr);
        vmr_raw{p}= vmr{p};
        vmr{p}    = humidity_transform(RH{p},'RH','VMR',T{p},pres{p});
        rho_raw{p}= rho{p};
        rho{p}    = humidity_transform(RH{p},'RH','HUM',T{p},pres{p});

        % calculate iwv
        ind     = find(~isnan(rho{p})); % humidity data is available only up to 10-12 km
        if length(ind)>1
            iwv(p)  = trapz(alt{p}(ind),rho{p}(ind)); 
        end
    end
    rsd.vmr_raw  = vmr_raw;
    rsd.RH_raw   = RH_raw;
    rsd.rho_raw  = rho_raw;
    rsd.iwv_raw  = iwv_raw;
end

rsd.vmr       = vmr;
rsd.RH        = RH;
rsd.rho       = rho;
rsd.iwv       = iwv;
rsd.rho_air   = rho_air;
rsd.wdir      = wdir;
rsd.wv        = wv;
