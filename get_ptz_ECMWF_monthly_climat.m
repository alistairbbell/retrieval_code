% function ptz = get_ptz_ECMWF_monthly_climat(start_time_str, stop_time_str, lat, lon, month_number)
%
% 2003-03-13 bd
% History:
% 2005-04-18 ah: Wechsel auf ECMWF ptz.
% 2005-09-06 ah: Schreibprozess wurde in uebergeordnetes skript ausgelagert.
% 2006-04-13 ah: Neu 91 ECMWF levels.
% 2011-07-04 ds: Vereinheitlichung mit anderen get_ptz_ Routinen
% 2012-07-04 ds
%
function ptz = get_ptz_ECMWF_monthly_climat(start_time_str, stop_time_str, lat, lon, month_number)

ptz = [];

% Find closest grid point
lat_grid = -90:1.125:90;
[tmpvar,lat_grid_ind] = min(abs(lat-lat_grid));
lat_grid_pt = lat_grid(lat_grid_ind);
lon_grid=-180:1.125:178.875;
[tmpvar,lon_grid_ind] = min(abs(lon-lon_grid));
lon_grid_pt = lon_grid(lon_grid_ind);

% Query data from MySQL
ecmwf_data = mysql(['select round(to_days(date)+1) as time, hour, p, T, z from atmospheres.ECMWF where date between "' datestr(floor(datenum(start_time_str,31)-6/24),29) '" and "' datestr(floor(datenum(stop_time_str,31)+6/24),29) '" and month(date) = ' num2str(month_number) ' and latitude = ' num2str(lat_grid_pt) ' and longitude = ' num2str(lon_grid_pt) ';'],'mat');

if isempty(ecmwf_data)
    disp('No ECMWF data found!')
    return
end

ecmwf_data = sortrows(ecmwf_data,[1 2 3]);

% Determine number of levels and times
n_time = length(unique(ecmwf_data(:,1)+ecmwf_data(:,2)/24));
n_levels = round(size(ecmwf_data,1)/n_time);

if n_time*n_levels~=size(ecmwf_data,1)
    disp('ERROR: Mix-up of ECMWF data...')
    return
end

% Make a matrix out of the data vectors:
time = reshape(ecmwf_data(:,1)+ecmwf_data(:,2)/24,n_levels,n_time);
p = reshape(ecmwf_data(:,3),n_levels,n_time);
T = reshape(ecmwf_data(:,4),n_levels,n_time);
z = reshape(ecmwf_data(:,5),n_levels,n_time);

% Interpolate to constant pressure levels
p_i = power(10,nanmean(log10(p),2));
T_i = NaN(size(T));
z_i = NaN(size(z));
for i = 1:size(T,2)
    T_i(:,i) = interp1(log10(p(:,i)),T(:,i),log10(p_i));
    z_i(:,i) = interp1(log10(p(:,i)),z(:,i),log10(p_i));
end

p = p_i;
T = nanmean(T_i,2);
z = nanmean(z_i,2);
ptz = flipud([p T z]);
ptz = ptz(~isnan(ptz(:,1)) & ~isnan(ptz(:,2)) & ~isnan(ptz(:,3)),:);

% get CIRA data for heights above 80 km
[T_cira, D_cira] = get_cira86(datenum(2010,month_number,1), lat);

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
ptz=[ptz; D_cira_i(ind) T_cira_i(ind) Z_cira_i(ind)];


% get USSTD data for the lowest heights
ptz_usstd_raw = [[101325 22632.1 5474.89 868.019 110.906 66.9389 3.95642 0.3734 0.01]',[288.15 216.65 216.65 228.65 270.65 270.65 214.65 273.15-86.2 273.15-86.2]',[0 11000 20000 32000 47000 51000 71000 84852 100000]'];
clear ptz_usstd
ptz_usstd(:,1) = power(10,linspace(log10(ptz_usstd_raw(1,1)),log10(ptz_usstd_raw(end,1)),100))';
ptz_usstd(:,2) = interp1(log10(ptz_usstd_raw(:,1)),ptz_usstd_raw(:,2),log10(ptz_usstd(:,1)));
ptz_usstd(:,3) = interp1(log10(ptz_usstd_raw(:,1)),ptz_usstd_raw(:,3),log10(ptz_usstd(:,1)));
ptz_usstd(:,2) = smooth(ptz_usstd(:,2),5);
ptz_usstd(:,3) = smooth(ptz_usstd(:,3),5);

% find USSTD indices below ECMWF profile
ind=find(ptz_usstd(:,3)<min(ptz(:,3)) & ptz_usstd(:,1)>max(ptz(:,1)));

% merge the two profiles
ptz=[ptz_usstd(ind,:); ptz];
