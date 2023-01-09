% function mlsdata = get_ptz_MLS_tempera_from_mysql(start_time_str, stop_time_str, lat, lon, mls_vers)
%
% 1. column = Pressure
% 2. column = Temperature partly from Tempera, where MR is high enough.
% 3. column = GPH
%
% Martin Lainer, 2015-07-21, Update to Aura MLS v4.2
% Martin Lainer, 2015-09-25, Tempera inclusion

function mlsdata = get_ptz_MLS_tempera_from_mysql(start_time_str, stop_time_str, lat, lon, mls_vers)

mlsdata = [];

species = 'Temperature';

if ~exist('mls_vers','var'); mls_vers = 4.2; end

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

% The time interval should be completely covered by MLS data and Tempera
% data

% get MLS and Tempera data
try mysql('drop table tmp.mls_profiles;'); end
mysql(['create temporary table tmp.mls_profiles select profile_ID,round(to_days(time)+1+time_to_sec(time)/86400,6) as time from satellites.EOS_MLS_header where time between "' datestr(start_time-4,31) '" and "' datestr(stop_time+4,31) '" and latitude between ' num2str(min(dlat)) ' and ' num2str(max(dlat)) ' and longitude between ' num2str(min(dlon)) ' and ' num2str(max(dlon)) ' and species = "GPH" and version = ' num2str(mls_vers) ' order by time;']);
GPH=mysql('select tmp.mls_profiles.time,pressure,value,satellites.EOS_MLS_profile.profile_ID from satellites.EOS_MLS_profile,tmp.mls_profiles where satellites.EOS_MLS_profile.profile_ID=tmp.mls_profiles.profile_ID','mat');
mysql('drop table tmp.mls_profiles;');

if isempty(GPH)
    disp('No Aura MLS data found!')
    return
end

mls_pressure = GPH(1:36,2); %GPH, T and H20 have the same pressure grid [36 x 1], only valid for MLS v4.2

% get Temperatur profile(s) from TEMPERA radiometer
Tempera = get_tempera_from_mysql(datestr(start_time,31),datestr(stop_time,31),4);

if isempty(Tempera)
    disp('No Tempera data found!')
    return
end

%Tempera.T_good = NaN(size(Tempera.T,1),size(Tempera.T,2));
%Tempera.p_good = NaN(size(Tempera.T,1),size(Tempera.T,2));

for i=1:size(Tempera.T,2)
    Tempera.T_good(:,i) = Tempera.T(Tempera.p(:,i) >= 90 & Tempera.p(:,i) <= 5000,i); %ca. 20 -50 km
    Tempera.p_good(:,i) = Tempera.p(Tempera.p(:,i) >= 90 & Tempera.p(:,i) <= 5000,i);
end

Tempera.T_mean = nanmean(Tempera.T_good,2);
Tempera.p = Tempera.p_good(:,1);

Tempera.T_mean_i = interp1(log(Tempera.p),Tempera.T_mean,log(mls_pressure));

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

VAR_ii_mean = nanmean(VAR_ii,2); % Same size as mls_pressure
VAR_TEMP_MLS = NaN(length(mls_pressure),1);
VAR_TEMP_MLS(~isnan(Tempera.T_mean_i)) = Tempera.T_mean_i(~isnan(Tempera.T_mean_i));
VAR_TEMP_MLS(isnan(Tempera.T_mean_i)) = VAR_ii_mean(isnan(Tempera.T_mean_i));

%Smooth Temperature profile
T_smooth = smooth(VAR_TEMP_MLS);

%if 0
%plot T profile
%figure('Visible','On');
%hold on;
%plot(VAR_TEMP_MLS,mls_pressure/100,'--r');
%plot(T_smooth,mls_pressure/100,'b');
%plot(nanmean(VAR_ii,2),mls_pressure/100,'--k')
%xlabel(gca,'T [K]');
%ylabel(gca,'Presuure [hPa]');
%title('2013-01-01 to 2013-01-02');
%legend('TEMPERA & MLS','TEMPERA & MLS (SMOOTH)','MLS')
%set(gca,'Yscale','log','YDir','reverse');
%end

mlsdata = [mls_pressure VAR_TEMP_MLS nanmean(GPH_ii,2)];
mlsdata = mlsdata(~isnan(mlsdata(:,2)) | ~isnan(mlsdata(:,3)),:);

mlsdata_smooth = [mls_pressure T_smooth nanmean(GPH_ii,2)];
mlsdata_smooth = mlsdata_smooth(~isnan(mlsdata_smooth(:,2)) | ~isnan(mlsdata_smooth(:,3)),:);

ptz = mlsdata;
ptz_smooth = mlsdata_smooth;

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
ind_smooth=find(Z_cira_i>max(ptz_smooth(:,3)) & Z_cira_i>95000 & D_cira_i<min(ptz_smooth(:,1)));

% merge the two profiles
mlsdata = [ptz; D_cira_i(ind) T_cira_i(ind) Z_cira_i(ind)];
%mlsdata_smooth = [ptz_smooth; D_cira_i(ind_smooth) T_cira_i(ind_smooth) Z_cira_i(ind_smooth)];
