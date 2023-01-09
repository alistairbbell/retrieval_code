% function xa=get_h2o_apr_MLSclimat(time,latitude)
%
% This function interpolates a zonal climatology to
% a given time and latitude. The zonal climatology
% is derive from MLS data and is used as a priori
% in the Aura/MLS retrieval for p<10 hPa.
% reference: Read et al. 2007, Section 2.3
%
% input: time as datenumber
%        latitude [-80 80]
%
% output: pressure and mixing ratio profile
%
%
function xa=get_h2o_apr_MLSclimat(time,latitude)

% get data from mysql
C=mysql('select day_of_year,latitude,p,H2O_VMR from atmospheres.MLS_Aura_H2O_climatology','mat');

% reshape data
doy=reshape(C(:,1),85,33,12);
lat=reshape(C(:,2),85,33,12);
p=reshape(C(:,3),85,33,12);
x=reshape(C(:,4),85,33,12);


% for proper interpolation append the first month at the end of data and the last month at the
% beginning: extended doy/data
edoy=nan(85,33,14);
edoy(:,:,1)=doy(:,:,1)-30;
edoy(:,:,2:13)=doy;
edoy(:,:,14)=doy(:,:,end)+30;

elat=nan(85,33,14);
elat(:,:,1)=lat(:,:,end);
elat(:,:,2:13)=lat;
elat(:,:,14)=lat(:,:,1);

ep=nan(85,33,14);
ep(:,:,1)=p(:,:,end);
ep(:,:,2:13)=p;
ep(:,:,14)=p(:,:,1);

ex=nan(85,33,14);
ex(:,:,1)=x(:,:,end);
ex(:,:,2:13)=x;
ex(:,:,14)=x(:,:,1);


% prepare interpolation (convert time to day-of-year)
p_i=p(:,1,1);
lat_i=latitude*ones(size(p_i));
doy_i=(time-datenum(year(time),1,0))*ones(size(p_i));

% prepare output
xa=nan(length(p_i),2);

% do interpolation
xa(:,1)=p_i;
xa(:,2)=interp3(elat,ep,edoy,ex,lat_i,p_i,doy_i);