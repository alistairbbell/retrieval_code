% get_avk_params(avk,MR,p,z)
%
% Check avks to get some parameters, which are inserted into mysql-DB:
%
% in:      avk             averaging kernel
%          MR              measurement response
%          p, z            p and z grid
%
% out:     ap.MRmax           maximal value of "measurement response"
%          ap.zMRmax          altitude of MRmax
%          ap.MR08b           bottom altitude of MRmax>=0.8
%          ap.MR08t           top altitude of MRmax>=0.8
%
%          ap.z_nom           nominal altitude
%          ap.z_avkmax        altitude of AVK max
%          ap.avkmax          AVK max value
%          ap.z_avkmaxf       altitude of AVK max (peak of 3th order
%                             polynomial fit)
%          ap.avkmaxf         AVK max value (fit)
%          ap.FWHM            Full width at half-maximum of each AVK
%
% (C) Rene Bleisch 10.2010, rene.bleisch@iap.unibe.ch
function ap = get_avk_params(avk,MR,p,z,nzz)
fid = 1;
nz  = length(p);

if nargin==4
    nzz  = 11;
    ntop = nzz+2;
else
    ntop = nzz;
end

% max MR / MR>0.8
[MRmax,ind] = nanmax(MR);
zMRmax      = z(ind);
     
MR08b       = interp1(MR(1:ind),z(1:ind),0.8);
MR08t       = interp1(MR(ind:end),z(ind:end),0.8);

% AVK max, FWHM
avkmax      = nan(nzz,1);
zavkmax     = nan(nzz,1);   
fwhm        = nan(nzz,1);
avkmaxf     = nan(nzz,1);
zavkmaxf    = nan(nzz,1); 
fwhmf       = nan(nzz,1);

% Disable warnings to omit lots of output 
warning off all;

for i = 1:nzz
    % Max
    [avkmax(i),ind0] = max(avk(i,:));
    zavkmax(i)      = z(ind0);
        
    % FWHM
    try
        zhmb     = interp1(avk(i,1:ind0),z(1:ind0),avkmax(i)/2); 
    catch
        eq       = find(diff(avk(i,1:ind0))==0);
        tmpi     = 1:ind0;
        fprintf(fid,'Problems during get_avk_params (L57), check AVK row %d / col %d/%d\n',i,tmpi(eq),tmpi(eq+1));
        continue
    end

    try
        zhmt     = interp1(avk(i,ind0:ntop),z(ind0:ntop),avkmax(i)/2); 
    catch
        eq       = find(diff(avk(i,ind0:ntop))==0);
        tmpi     = ind0:nzz+2;
        fprintf(fid,'Problems during get_avk_params (L66), check AVK row %d / col %d/%d\n',i,tmpi(eq),tmpi(eq+1));
        continue
    end
    fwhm(i)  = zhmt-zhmb;
        
    % Fit a 3rd order polynom
    ind     = find(avk(i,:)>1/2*avkmax(i));
    if length(ind)<3
        if ind0>2
            ind = ind0-2:ind0+2;
        elseif ind0>1
            ind = ind0-1:ind0+2;
        else
            ind = ind0:ind0+2;
        end
    end
    pf      = polyfit(z(ind)',avk(i,ind),3);
    zpf     = z(ind(1)):10:z(ind(end));
    avkp    = polyval(pf,zpf);
    [avkmaxf(i),ind] = max(avkp);
    zavkmaxf(i) = zpf(ind);
        
    if avkmaxf(i)<avkmax(i)
        avkmaxf(i) = avkmax(i);
        zavkmaxf(i) = zavkmax(i);
    end
        
    % FWHM fit
    zhmbf     = interp1(avk(i,1:ind0),z(1:ind0),avkmaxf(i)/2); 
    zhmtf     = interp1(avk(i,ind0:ntop),z(ind0:ntop),avkmaxf(i)/2); 

    fwhmf(i)  = zhmtf-zhmbf;        
end

% Enable warnings!
warning on all;

%1D-data
ap.MRmax   = MRmax;
ap.zMRmax  = zMRmax;
ap.MR08b   = MR08b;
ap.MR08t   = MR08t;

%2D-data
ap.z_nom      = z(1:nzz);
ap.z_avkmax   = zavkmax;
ap.avkmax     = avkmax;
ap.z_avkmaxf  = zavkmaxf;
ap.avkmaxf    = avkmaxf;
ap.fwhm       = fwhmf;


end