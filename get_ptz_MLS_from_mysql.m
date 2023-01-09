% function mlsdata = get_ptz_MLS_from_mysql(start_time_str, stop_time_str, lat, lon, mls_vers)
%
% 1. column = Pressure
% 2. column = Temperature
% 3. column = GPH
%
% Dominik Scheiben, 2011-07-04
% Martin Lainer, 2015-07-21, Update to Aura MLS v4.2

function mlsdata = get_ptz_MLS_from_mysql(start_time_str, stop_time_str, lat, lon, mls_vers)

% start_time_str = '2011-01-01';
% stop_time_str = '2011-01-10';
% lat = 47;
% lon = 7;
% mls_vers = 3.3;

mlsdata = [];

species = 'Temperature';

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

if mls_vers == 4.2
    mls_pressure = power(10,[(-1:1/3:2/3)'; (1:1/6:11/6)'; (2:1/12:12/3)']); %Consider LPD change per decade in pressure, according to MLS data quality and description doc.
    mls_pressure = mls_pressure(1:end-1);
elseif mls_vers == 3.3
    mls_pressure = power(10,[(-1:1/3:2/3)'; (1:1/6:11/6)'; (2:1/12:12/3)'; 12.9/3]);
elseif mls_vers == 2.2
    mls_pressure = power(10,[(-2/3:1/3:1)'; (7/6:1/6:9.5/3)'; (20/6:1/12:26.0/6)'; 26.4/6]);
end


% The time interval should be completely covered by MLS data


% get MLS data
try mysql('drop table tmp.mls_profiles;'); end
mysql(['create temporary table tmp.mls_profiles select profile_ID,round(to_days(time)+1+time_to_sec(time)/86400,6) as time from satellites.EOS_MLS_header where time between "' datestr(start_time-4,31) '" and "' datestr(stop_time+4,31) '" and latitude between ' num2str(min(dlat)) ' and ' num2str(max(dlat)) ' and longitude between ' num2str(min(dlon)) ' and ' num2str(max(dlon)) ' and species = "GPH" and version = ' num2str(mls_vers) ' order by time;']);
GPH=mysql('select tmp.mls_profiles.time,pressure,value,satellites.EOS_MLS_profile.profile_ID from satellites.EOS_MLS_profile,tmp.mls_profiles where satellites.EOS_MLS_profile.profile_ID=tmp.mls_profiles.profile_ID','mat');
mysql('drop table tmp.mls_profiles;');
if isempty(GPH)
    disp('No Aura MLS data found!')
    return
end
mls_times_GPH_i = floor(min(GPH(:,1))*2)/2+0.25:0.5:ceil(max(GPH(:,1))*2)/2-0.25;
GPH_i = NaN(length(mls_pressure),length(mls_times_GPH_i));
for i = 1:length(mls_times_GPH_i)
    profile_IDs = unique(GPH(GPH(:,1)>=mls_times_GPH_i(i)-0.25 & GPH(:,1)<mls_times_GPH_i(i)+0.25,4));
    if ~isempty(profile_IDs)
        value_i = NaN(length(mls_pressure),length(profile_IDs));
        for j = 1:length(profile_IDs)
            ind = find(GPH(:,4)==profile_IDs(j));
            if ~isempty(ind)
                pressure_tmp = GPH(ind,2);
                value_tmp = GPH(ind,3);
                value_i(:,j) = interp1(log(pressure_tmp),value_tmp,log(mls_pressure));
            end
        end
        GPH_i(:,i) = nanmean(value_i,2);
    end
end
mls_times_GPH_i = mls_times_GPH_i(~isnan(nanmean(GPH_i,1)));
GPH_i = GPH_i(:,~isnan(nanmean(GPH_i,1)));
if max(mls_times_GPH_i)<stop_time || min(mls_times_GPH_i)>stop_time
    disp('Aura MLS data does not cover the whole requested time period!');
    return
end
mysql(['create temporary table tmp.mls_profiles select profile_ID,round(to_days(time)+1+time_to_sec(time)/86400,6) as time from satellites.EOS_MLS_header where time between "' datestr(start_time-4,31) '" and "' datestr(stop_time+4,31) '" and latitude between ' num2str(min(dlat)) ' and ' num2str(max(dlat)) ' and longitude between ' num2str(min(dlon)) ' and ' num2str(max(dlon)) ' and species = "' species '" and version = ' num2str(mls_vers) ' order by time;']);
VAR=mysql('select tmp.mls_profiles.time,pressure,value,satellites.EOS_MLS_profile.profile_ID from satellites.EOS_MLS_profile,tmp.mls_profiles where satellites.EOS_MLS_profile.profile_ID=tmp.mls_profiles.profile_ID','mat');
mysql('drop table tmp.mls_profiles;');
if isempty(VAR)
    disp('No Aura MLS data found!')
    return
end
mls_times_VAR_i = floor(min(VAR(:,1))*2)/2+0.25:0.5:ceil(max(VAR(:,1))*2)/2-0.25;
VAR_i = NaN(length(mls_pressure),length(mls_times_VAR_i));
for i = 1:length(mls_times_VAR_i)
    profile_IDs = unique(VAR(VAR(:,1)>=mls_times_VAR_i(i)-0.25 & VAR(:,1)<mls_times_VAR_i(i)+0.25,4));
    if ~isempty(profile_IDs)
        value_i = NaN(length(mls_pressure),length(profile_IDs));
        for j = 1:length(profile_IDs)
            ind = find(VAR(:,4)==profile_IDs(j));
            if ~isempty(ind)
                pressure_tmp = VAR(ind,2);
                value_tmp = VAR(ind,3);
                value_i(:,j) = interp1(log(pressure_tmp),value_tmp,log(mls_pressure));
            end
        end
        VAR_i(:,i) = nanmean(value_i,2);
    end
end
mls_times_VAR_i = mls_times_VAR_i(~isnan(nanmean(VAR_i,1)));
VAR_i = VAR_i(:,~isnan(nanmean(VAR_i,1)));
if max(mls_times_VAR_i)<stop_time || min(mls_times_VAR_i)>stop_time
    disp('ERROR: Aura MLS data does not cover the whole requested time period!');
    return
end

times_mls_ii = start_time:1/24:stop_time;
GPH_ii = interp2(mls_times_GPH_i,log(mls_pressure),GPH_i,times_mls_ii,log(mls_pressure));
VAR_ii = interp2(mls_times_VAR_i,log(mls_pressure),VAR_i,times_mls_ii,log(mls_pressure));

mlsdata = [mls_pressure nanmean(VAR_ii,2) nanmean(GPH_ii,2)];
mlsdata = flipud(mlsdata(~isnan(mlsdata(:,2)) | ~isnan(mlsdata(:,3)),:));

ptz = mlsdata;

% get CIRA data for heights above 80 km
[T_cira, D_cira] = get_cira86(mean([datenum(start_time_str,31) datenum(stop_time_str,31)]), lat);

Z_cira_i=(0:1000:120000)';
T_cira_i=interp1(T_cira(:,1),T_cira(:,2),Z_cira_i);
D_cira_i=exp(interp1(D_cira(:,1),log(D_cira(:,2)),Z_cira_i));

ind=find(isnan(T_cira_i)==0 & isnan(D_cira_i)==0);

Z_cira_i=Z_cira_i(ind);
T_cira_i=T_cira_i(ind);
D_cira_i=D_cira_i(ind);

% find CIRA indices above ECMWF profile
% ind=find(Z_cira_i>ptz(end,3) & D_cira_i<ptz(end,1));
ind=find(Z_cira_i>max(ptz(:,3)) & Z_cira_i>95000 & D_cira_i<min(ptz(:,1)));

% merge the two profiles
mlsdata = [ptz; D_cira_i(ind) T_cira_i(ind) Z_cira_i(ind)];


