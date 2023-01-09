% function ptz = get_ptz_ECMWF_temporal_interpolation(start_time_str, stop_time_str, lat, lon)
%
% 2003-03-13 bd
% History:
% 2005-04-18 ah: Wechsel auf ECMWF ptz.
% 2005-09-06 ah: Schreibprozess wurde in uebergeordnetes skript ausgelagert.
% 2006-04-13 ah: Neu 91 ECMWF levels.
% 2011-07-04 ds: Vereinheitlichung mit anderen get_ptz_ Routinen
%
function ptz = get_ptz_ECMWF_new(start_time_str, stop_time_str, lat, lon)

ptz = [];


YYYY = datestr(start_time_str,'yyyy');
MM   = datestr(start_time_str,'mm');
DD   = datestr(start_time_str,'dd');

file_ecmwf_nya = ['/export/data/miawarac/ECMWF_DATA_NYA/ECMWF_OPER_v2_' YYYY MM DD '_NYA_cat.nc'];

if ~exist(file_ecmwf_nya,'file')
    file_ecmwf = [ '/storage/tub/atmosphere/ecmwf/oper/' YYYY '/ECMWF_OPER_v2_' YYYY MM DD '.nc'];
    date_num = floor(datenum(start_time_str));
    extract_specific_location_from_ecmwf_nc(date_num,lat,lon,file_ecmwf,file_ecmwf_nya)
end



time = (double(squeeze(ncread(file_ecmwf_nya,'time')))'/24)+datenum('1900-01-01'); % time in matlab format
lnsp = squeeze(ncread(file_ecmwf_nya,'lnsp')); % log of surface pressure 
T    = squeeze(ncread(file_ecmwf_nya,'t')); % Temperature
phi  = squeeze(ncread(file_ecmwf_nya,'z')); % surface geopotential
q    = squeeze(ncread(file_ecmwf_nya,'q')); % specific humidity

%addpath('/home/gromosc/matlab/GromosC/datafiles')
load('ecmwf_hybrid_levels_137.mat')

p_surf = exp(lnsp(1,:));
Phi_surf = phi(1,:);
for k = 1:length(time)
	[ z(:,k), p(:,k) ] = ecmwf_zp_calc( A_h, B_h, T(:,k), p_surf(k), Phi_surf(k), q(:,k) );
end

% for compatibility with database use same data format
z = round(z);
p = round(p);
T = round(T);

% Interpolate to constant pressure levels
p_i = power(10,nanmean(log10(p),2));
T_i = NaN(size(T));
z_i = NaN(size(z));
for i = 1:size(T,2)
    T_i(:,i) = interp1(log10(p(:,i)),T(:,i),log10(p_i));
    z_i(:,i) = interp1(log10(p(:,i)),z(:,i),log10(p_i));
end

% Interpolate to 1-hourly time vector and then take the average
time_i = datenum(datestr(datenum(start_time_str,31),'yyyymmddHH'),'yyyymmddHH'):1/24:datenum(datestr(datenum(stop_time_str,31)+1/24-3/3600/24,'yyyymmddHH'),'yyyymmddHH');
T_i2 = interp2(time(1,:),log(p_i),T_i,time_i,log(p_i));
z_i2 = interp2(time(1,:),log(p_i),z_i,time_i,log(p_i));
p = p_i;
T = nanmean(T_i2,2);
z = nanmean(z_i2,2);
ptz = flipud([p T z]);
ptz = ptz(~isnan(ptz(:,1)) & ~isnan(ptz(:,2)) & ~isnan(ptz(:,3)),:);

if isempty(ptz)
    return
end

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
k
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
