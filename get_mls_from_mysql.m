% function mlsdata = get_mls_from_mysql(start_time_str, stop_time_str, lat, lon, species, mls_vers)
%
% 1. column = Pressure
% 2. column = Species
% 3. column = GPH
%
%
% Dominik Scheiben, 2011-07-04
%
function mlsdata = get_mls_from_mysql(start_time_str, stop_time_str, lat, lon, species, mls_vers)

% start_time_str = '2011-01-01';
% stop_time_str = '2011-01-10';
% lat = 47;
% lon = 7;
% species = 'H2O';
% mls_vers = 3.3;

mlsdata = [];

if ~exist('mls_vers','var'); mls_vers = 3.3; end

start_time = datenum(start_time_str,31);
stop_time = datenum(stop_time_str,31);

% meridional, zonal tolerance [m]
[dlat, dlon] = satellite_area(lat,lon,400,800);

% the time interval should be 3 days or more
if stop_time-start_time<3
    t1 = mean([stop_time start_time])-1.5;
    t2 = mean([stop_time start_time])+1.5;
    start_time = t1;
    stop_time  = t2;
end

% get MLS data
try mysql('drop table tmp.mls_profiles;'); end
mysql(['create temporary table tmp.mls_profiles select profile_ID,round(to_days(time)+1+time_to_sec(time)/86400,6) as time from satellites.EOS_MLS_header where time between "' datestr(start_time,31) '" and "' datestr(stop_time,31) '" and latitude between ' num2str(min(dlat)) ' and ' num2str(max(dlat)) ' and longitude between ' num2str(min(dlon)) ' and ' num2str(max(dlon)) ' and species = "GPH" and version = ' num2str(mls_vers) ' order by time;']);
GPH=mysql('select tmp.mls_profiles.time,pressure,value,satellites.EOS_MLS_profile.profile_ID from satellites.EOS_MLS_profile,tmp.mls_profiles where satellites.EOS_MLS_profile.profile_ID=tmp.mls_profiles.profile_ID','mat');
mysql('drop table tmp.mls_profiles;');
mysql(['create temporary table tmp.mls_profiles select profile_ID,round(to_days(time)+1+time_to_sec(time)/86400,6) as time from satellites.EOS_MLS_header where time between "' datestr(start_time,31) '" and "' datestr(stop_time,31) '" and latitude between ' num2str(min(dlat)) ' and ' num2str(max(dlat)) ' and longitude between ' num2str(min(dlon)) ' and ' num2str(max(dlon)) ' and species = "' species '" and version = ' num2str(mls_vers) ' order by time;']);
VAR=mysql('select tmp.mls_profiles.time,pressure,value,satellites.EOS_MLS_profile.profile_ID from satellites.EOS_MLS_profile,tmp.mls_profiles where satellites.EOS_MLS_profile.profile_ID=tmp.mls_profiles.profile_ID','mat');
mysql('drop table tmp.mls_profiles;');

if isempty(GPH) || isempty(VAR)
    disp('No Aura MLS data found!')
    return
end

if mls_vers == 3.3
    mls_pressure = power(10,[(-1:1/3:2/3)'; (1:1/6:11/6)'; (2:1/12:13/3)']);
elseif mls_vers == 2.2
    mls_pressure = power(10,[(-2/3:1/3:1)'; (7/6:1/6:9.5/3)'; (20/6:1/12:26.5/6)']);
end

GPH = sortrows(GPH,[1 2]);
VAR   = sortrows(VAR,  [1 2]);

profiles_GPH = unique(GPH(:,4));
profiles_VAR = unique(VAR(:,4));

n_time_GPH = length(profiles_GPH);
n_time_VAR = length(profiles_VAR);

Z_i = NaN(length(mls_pressure),n_time_GPH);
for i = 1:length(profiles_GPH)
   Z_i(:,i) = interp1(log10(GPH(GPH(:,4)==profiles_GPH(i),2)),GPH(GPH(:,4)==profiles_GPH(i),3),log10(mls_pressure));
end
VAR_i = NaN(length(mls_pressure),n_time_VAR);
for i = 1:length(profiles_VAR)
   VAR_i(:,i) = interp1(log10(VAR(VAR(:,4)==profiles_VAR(i),2)),VAR(VAR(:,4)==profiles_VAR(i),3),log10(mls_pressure));
end

mlsdata = [mls_pressure nanmean(VAR_i,2) nanmean(Z_i,2)];
mlsdata = mlsdata(~isnan(mlsdata(:,2)),:);


