function [msm]=get_level1_acq_from_NetCDF_certain_time_period(start_date_str,end_date_str,pol0)



% filename

fileL1 = path_names(start_date_str);

if exist(fileL1,'file') 

	% ncread
	time             = ncread(fileL1,'/spectrometer1/time');
    max_time         = ncread(fileL1,'/spectrometer1/last_sky_time');
    min_time         = ncread(fileL1,'/spectrometer1/first_sky_time');
    altitude         = ncread(fileL1,'/spectrometer1/alt');
    latitude         = ncread(fileL1,'/spectrometer1/lat');
    longitude        = ncread(fileL1,'/spectrometer1/lon');
    tint             = ncread(fileL1,'/spectrometer1/calibration_time');
    f                = ncread(fileL1,'/spectrometer1/frequencies');
    y                = ncread(fileL1,'/spectrometer1/Tb');
    sigma            = ncread(fileL1,'/spectrometer1/sigma');
    %integration_ok   = ncread(fileL1,'/spectrometer1/time');
    tau              = ncread(fileL1,'/spectrometer1/mean_opacity');
    trec             = ncread(fileL1,'/spectrometer1/TSys');
    a                = ncread(fileL1,'/spectrometer1/A');
    A                = ncread(fileL1,'/spectrometer1/a');
    mirror_elevation = ncread(fileL1,'/spectrometer1/mean_sky_elevation_angle');
	pol          = ncread(fileL1,'/spectrometer1/pol');
    % select date & pol
    
    idx_t = time+datenum('1970-01-01') >= datenum(start_date_str) & time+datenum('1970-01-01') <= datenum(end_date_str);
    idx_pol =  pol == pol0;
    idx = idx_t & idx_pol;

    % output format
    msm.time             = mean(time(idx)+datenum('1970-01-01'));
    msm.max_time         = mean(max_time(idx)+datenum('1970-01-01'));
    msm.min_time         = mean(min_time(idx)+datenum('1970-01-01'));
    msm.altitude         = altitude(1);
    msm.latitude         = latitude(1);
    msm.longitude        = longitude(1);
    msm.tint             = tint(1);
    msm.f                = f;
    msm.y                = mean(y(:,idx),2)';
    %msm.sigma            = mean(sigma(idx));
    msm.integration_ok   = 1;
    msm.tau              = nanmean(tau(idx));
    msm.trec             = mean(trec(idx));
    msm.a                = mean(a(idx));
    msm.A                = mean(A(idx));
    msm.mirror_elevation = mean(mirror_elevation(idx));
    msm.pol		 = pol0;
   
    ind = find(f>=2.226002048709027e+10 & f>=2.226496463651346e+10);

    m = length(ind);
    [p,s,mu]=polyfit(1:m,msm.y(ind),1);
    level1_sigma=msm.y(ind)-polyval(p,1:m,[],mu);
    msm.sigma=std(level1_sigma);
    
    
else
	msm = empty_output;
end



%== subroutines ==%

function msm=empty_output
msm.time=[];
msm.max_time=[];
msm.min_time=[];
msm.altitude=[];
msm.latitude=[];
msm.longitude=[];
msm.tint=[];
msm.f=[];
msm.y=[];
msm.sigma=[];
msm.integration_ok=0;
msm.tint=[];
msm.tau=[];
msm.trec=[];
msm.a=[];
msm.A=[];
msm.mirror_elevation=[];
msm.pol= [];
