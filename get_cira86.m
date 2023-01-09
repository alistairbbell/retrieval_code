function [temp_profil, druck_profil] = get_cira86(time, latitude)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%    G  E  T  _  C  I  R  A  8  6  .  M                         %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
% Funktion, die aus der Cira86 - Mat.Datei das Druck- und       %
% Temperaturprofil zu einer verlangten Breite und Monat         %
% herauszieht.                                                  %
%                                                               %
% Input:        month         --  Monat  (MM)                   %
%               latitude      --  Breite (wird auf 10 Grad gerundet)%
%                                                               %
% Output:       temp_profil(:;1)  Hoehenvektor  [m]              %
%               temp_profil(:,2)  Temperaturen [K]              %
%               druck_profil(:,1) Hoehenvektor  [m]              %
%               druck_profil(:,2) Druck        [Pa]             %
%                                                               %
% A. Siegenthaler, September 1999                               %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> Processing cira86 data');

time_vec = datevec(time);

month = time_vec( 2 );

load cira86.mat;
% zuordnen der Breite auf einen der Werte:
% [-80 -70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 70 80]

breite = (round(latitude / 10 )) * 10;
if breite > 80
  breite = 80;
end
if breite < -80
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



