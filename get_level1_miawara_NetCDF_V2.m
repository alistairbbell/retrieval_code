function [msm]=get_level1_miawara_NetCDF_V2(time_str,extra_info)

%==========================================================================
% NAME      | get_level1_miawara_NetCDF_V2.m
% TYPE      | Script
% AUTHOR(S) | Alistair Bell
% CREATION  | 03.2023
%           |
% ABSTRACT  | Function to read in the measurment data from a NetCDF file.
%           | Data is then output to the msm structure 
%           | 
% ARGUMENTS | INPUTS: - time_str
%           |         - extra_info (structure)
%           |
%           | OUTPUTS: msm
%           |
%           |
%==========================================================================

% filename
fileL1 = join([extra_info.l1_dir, extra_info.l1_savename , time_str , extra_info.l1_file_ext]) ;
disp(fileL1)
pol0 = 0;

fc = extra_info.f_cen;
bw = extra_info.bw_obs;

if exist(fileL1,'file')
    disp('Found level 1 file')
	%ncread
	time             = ncread(fileL1,'/spectrometer1/time');
    max_time         = ncread(fileL1,'/spectrometer1/last_sky_time');
    min_time         = ncread(fileL1,'/spectrometer1/first_sky_time');
    altitude         = ncread(fileL1,'/spectrometer1/alt');
    latitude         = ncread(fileL1,'/spectrometer1/lat');
    longitude        = ncread(fileL1,'/spectrometer1/lon');
    tint             = ncread(fileL1,'/spectrometer1/calibration_time');
    f                = ncread(fileL1,'/spectrometer1/frequency');
    tau              = ncread(fileL1,'/spectrometer1/tau');
    y                = ncread(fileL1,'/spectrometer1/Tb'); 
    trec             = ncread(fileL1,'/spectrometer1/TSys');
    a                = ncread(fileL1,'/spectrometer1/A');
    A                = ncread(fileL1,'/spectrometer1/a');
    mirror_elevation = ncread(fileL1,'/spectrometer1/mean_sky_elevation_angle');
    meteo_time = ncread(fileL1,'/meteo/time');
    pressure = ncread(fileL1,'/meteo/air_pressure');

    % Convert UNIX time to MATLAB datenum
    unix_time = seconds(time);
    ref_date = datenum('1970-01-01');
    time_datenum = datenum(unix_time + ref_date);
    
    % Convert time_str to datenum
    time_str_datenum = datenum(time_str);
    
    % Create a 24-hour range for comparison
    time_range_start = time_str_datenum;
    time_range_end = time_str_datenum + hours(24);
    
    % Find the indices within the time range
    idx = (time_datenum >= time_range_start) & (time_datenum <= time_range_end);
    
    disp('y size')
    disp(size(y))
    
    disp('size(longitude)')
    disp(size(longitude))
    disp(size(latitude))
    disp(size(altitude))

 
    % output format
    for i = 1:length(time)
        msm(i).time             = time(i)/(3600*24);%mean(time(idx)+datenum('1970-01-01'));
        msm(i).max_time         = max_time(i)/(3600*24);%mean(max_time(idx)+datenum('1970-01-01'));
        msm(i).min_time         = min_time(i)/(3600*24);%mean(min_time(idx)+datenum('1970-01-01'));
        msm(i).altitude         = altitude(i);
        msm(i).latitude         = latitude(i);
        msm(i).longitude        = longitude(i);
        msm(i).tint             = tint(i);
        msm(i).f                = f(:);
        msm(i).y                = y(:,i);
        
        msm(i).integration_ok   = 1;
        msm(i).tau              = tau(i);
        msm(i).trec             = trec(i);
        msm(i).a                = mean(a(:));
        msm(i).A                = mean(A(:));
        msm(i).mirror_elevation = mean(mirror_elevation(:));
        msm(i).pol		 = pol0;
        msm(i).HSE_P = pressure(:);
        msm(i).meteo_time = meteo_time(:)/(3600*24);
 
        ind = find((f>=fc-bw/2) & (f>=fc-bw/2));
    
        m = length(ind);
         %[p,s,mu]=polyfit(1:m,msm(i).y(ind),1);
         %level1_sigma=msm(i).y(ind)'-polyval(p,1:m,[],mu);
         %disp('size(level1_sigma)')
         %disp(size(level1_sigma))
        end
else
    disp('File does not exist');
	msm = empty_output;
end

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
msm.HSE_P = [];
end