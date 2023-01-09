function ptz=get_pT_AuraMLS_from_mysql(time1, time2, lat, lon)
% function get_pT_ecmwf(time, lat, lon, z);
%
% Erstelle ptz Matrix
%
% time1:   Start Zeit als datenum
% time2:   End Zeit als datenum
% lat:     geo Breite
% lon:     geo Länge
%
% 2003-03-13 bd
% History:
% 2005-04-18 ah: Wechsel auf ECMWF ptz.
% 2005-09-06 ah: Schreibprozess wurde in übergeordnetes skript ausgelagert.
% 2006-04-13 ah: Neu 90 ECMWF levels.
% 2008-07-14 ah: Umgeschrieben für AuraMLS Temperatur Daten

% some constants
RE=6378e3; % [m] Earth radius

% meridional, zonal tolerance [m]
dx=400e3;
dy=800e3;
dlat=lat+[-1 1]*dx/RE/2/pi*360;
dlon=lon+[-1 1]*dy/(RE*cos(lat/360*2*pi))/2/pi*360;

% get MLS data
A=get_mls_from_mysql(datestr(time1,31),datestr(time2,31),dlat,dlon,'H2O',2.2);

% return, if there is no data
if isempty(A)==1
    disp('no MLS temperature profiles found');
    ptz=[];
    return
end

ptz(:,1)=A.p;
ptz(:,2)=nanmean(A.T')';
ptz(:,3)=nanmean(A.GPH')';


% remove nan's
ind=find(isnan(ptz(:,2))==0);
ptz=ptz(ind,:);
ind=find(isnan(ptz(:,3))==0);
ptz=ptz(ind,:);

if isempty(ptz)
    disp('no Aura/MLS data found')
    ptz=[];
    return
end

% get CIRA data
[T_cira, D_cira] = get_cira86(mean([time1 time2]), lat);

Z_cira_i=[0:1000:120000]';
T_cira_i=interp1(T_cira(:,1),T_cira(:,2),Z_cira_i);
D_cira_i=exp(interp1(D_cira(:,1),log(D_cira(:,2)),Z_cira_i));

ind=find(isnan(T_cira_i)==0 & isnan(D_cira_i)==0);

Z_cira_i=Z_cira_i(ind);
T_cira_i=T_cira_i(ind);
D_cira_i=D_cira_i(ind);

% find CIRA indices above AuraMLS profile
ind=find(Z_cira_i>ptz(end,3) & D_cira_i<ptz(end,1));

% merge the two profiles
ptz=[ptz; D_cira_i(ind) T_cira_i(ind) Z_cira_i(ind)];




%%%%%%%%%%%%%%%%%%% SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [temp_profil, druck_profil] = get_cira86(time, latitude)

disp('->     Processing cira86 data ...');

time_vec = datevec(time);

month = time_vec( 2 );

load /data/iapmw/atmos/database/atmospheres/cira.86/cira86.mat;

% zuordnen der Breite auf einen der Werte:
% [-80 -70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 70 80]

breite = (round(latitude / 10 )) * 10;
if breite > 80,
  breite = 80;
end
if breite < -80,
  breite = -80;
end

Temp = [cira86.T];
Druck = [cira86.p];
T_alt = [cira86.T_alt]*1000;
p_alt = [cira86.p_alt]*1000;

lat = [ -80 -70 -60 -50 -40 -30 -20 -10 0 ...
       10 20 30 40 50 60 70 80 ];
lat_index = find( breite == lat );

profil_index = ( month - 1 ) * length( lat ) + lat_index;

temp_profil(:,2) = Temp(:,profil_index);
temp_profil(:,1) = T_alt(:,month);

druck_profil(:,2) = Druck(:,profil_index);
druck_profil(:,1) = p_alt(:,month);
