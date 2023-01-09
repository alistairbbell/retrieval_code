function [ z, p ] = ecmwf_zp_calc( A_h, B_h, T, p_surf, Phi_surf, q )

% [ z, p ] = ecmwf_zp_calc( A_h, B_h, T, p_surf, Phi_surf, q )
%
% Calculates geometric altitude and pressure from ECMWF profiles on model
% levels.
%
% Input parameters:
%
% A_h:      A_k+1/2 constant from GDS section of GRIB file on half levels
% B_h:      B_k+1/2 constant from GDS section of GRIB file on half levels
% T:        Temperature [K] on full levels
% p_surf:   Surface pressure [Pa] (optional)
% Phi_surf: Surface geopotential [m^2 s^-2] (optional)
% q:        Specific humidity [kg kg^-1] on full levels (optional)
%
% If missing, the parameters p_surf, Phi_surf, and q will be initialized
% to p_surf = 1.01325e5 Pa, Phi_surf = 0 m^2/s^2, q = 0 kg/kg.
%
% Output parameters:
%
% z: geometric altitude [m] on full levels
% p: pressure [Pa] on full levels
%
% Note: the function assumes that all vectors are sorted so that the first
%       value is near the top of the atmosphere and the last value is near
%       the surface.
%
% $Id: ecmwf_zp_calc.m,v 1.1 2005/03/04 23:06:53 feist Exp $

% Nomenclature for full and half-level variables:
%
% X is a variable on full levels (like all the prognostic variables)
% X_h is a variable on half levels
%
% In the ECMWF documentation in chapter 2.2 that would be written as
% X   <-> X_k
% X_h <-> X_k+1/2 (which is not possible in Matlab syntax)

% Set missing variables to default values
if ~exist('p_surf', 'var')
  p_surf = 100*1013.250; % Pa
end

if ~exist('Phi_surf', 'var')
  Phi_surf = 0;
end

if ~exist('q', 'var')
  q = zeros(size(T)); % Specific humidity [kg/kg]
end

%
% Define physical constants according to the subroutine SUCST.F90
% from the source code of the ECMWF ITS model
%
RKBOL = 1.380658E-23;	% Boltzmann's constant k [J/K]
RNAVO = 6.0221367E+23;	% Avogadro's number NA []
R = RNAVO * RKBOL;	% ideal gas constant R [J/(mol*K)]
RMD = 28.9644;		% Dry air molecular weight [g/mol]
RMV = 18.0153;          % Water vapor molecular weight [g/mol]
RD = 1000 * R / RMD;	% Dry air constant Rd [J/(K*kg)]
RV = 1000 * R / RMV;    % Water vapor constand Rv [J/(K*kg)]
RG = 9.80665;		% Earth's gravitational acceleration g [m/s^2]

% Get number of levels
NLEV = length(T);

%
% Calculate delta_p according to Eq. 2.13
%
% Calculate half-level pressure p_h  from
% A_h, B_h and surface pressure p_surf
% Level 1 = top, NLEV = bottom
%
p_h = A_h + B_h * p_surf;              % Half-level pressure (eq. 2.11)


% Calculate delta_p according to Eq. 2.13

delta_p = diff( p_h );

%
% Calculate virtual temperature according to standard textbook formula
%
T_v = T .* ( 1 + ( RV / RD - 1 ) * q ); % Virtual temperature [K]

%
% Calculate ln_p = log( p_k+1/2 / p_k-1/2 ) for 1 <= k <= NLEV
%
S = warning; warning('off'); % Ignore warnings for this part
ln_p = log( p_h(2:end) ./ p_h(1:end-1) );
ln_p(1) = 0; % ln_p(1) is never used and might be infinite
warning(S);  % Reset previous warning state

%
% Geopotential height at half levels (eq. 2.21 in matrix form)
%
loopj = triu( ones( NLEV+1, NLEV ));  % Loop over j=k+1...NLEV as a matrix
Phi_h = Phi_surf + RD * loopj * ( T_v .* ln_p );
Phi_h(1) = NaN; % Phi_h(1) was calculated with ln_p(1) and is never used

%
% Calculate alpha according to eq. 2.23
%
alpha = 1 - p_h(1:end-1) ./ delta_p .* ln_p;
alpha(1) = log(2); % Top level specially defined

%
% Calculate geopotential on full levels according to eq. 2.22
%
Phi = Phi_h(2:end) + RD * alpha .* T_v;

%
% Return geometric height and pressure on full levels
%
p = ( p_h(1:end-1) + p_h(2:end) ) / 2; % Full level pressure [Pa]
z = Phi / RG; % Geometric height [m]

