% Get pTz-data from ECMWF and merge with surface data from ZIMM meteostation
% vmr is also included, as it is sometimes also used (simulations)
function [ptz,vmr] = get_pTz_ECMWF_trop(starttime,stoptime,lat,lon,do_cira)
if nargin==4
    do_cira = 0;
end

% find closest grid point
lat_grid = -90:1.125:90;
[tmpvar,lat_grid_ind] = min(abs(lat-lat_grid));
lat_grid_pt = lat_grid(lat_grid_ind);
lon_grid=-180:1.125:178.875;
[tmpvar,lon_grid_ind] = min(abs(lon-lon_grid));
lon_grid_pt = lon_grid(lon_grid_ind);

% Query data from MySQL
ecmwf_data = mysql(['select round(to_days(date)+1) as time, hour, p, T, z, H2O_vmr from atmospheres.ECMWF where date between "'...
    datestr(floor(starttime),29) '" and "' datestr(floor(stoptime),29) ...
    '" and latitude = ' num2str(lat_grid_pt) ' and longitude = ' num2str(lon_grid_pt) ';'],'mat');

if isempty(ecmwf_data)
    disp('No ECMWF data found!')
    ptz = [];
    vmr = [];
    return
end

ecmwf_data = sortrows(ecmwf_data,[1 2 3]);

time = unique(ecmwf_data(:,1)+ecmwf_data(:,2)/24);

time_ind_91 = find(time<datenum(2013,6,25,3,0,0));
ecmwf_data_91 = ecmwf_data(1:length(time_ind_91)*91,:);
time_ind_137 = find(time>datenum(2013,6,25,3,0,0));
ecmwf_data_137 = ecmwf_data(length(time_ind_91)*91+1:length(time_ind_91)*91+length(time_ind_137)*137,:);

if ~isempty(ecmwf_data_91) && ~isempty(ecmwf_data_137)
   
    % Determine number of levels and times with 91 levels:
    n_time_91 = length(unique(ecmwf_data_91(:,1)+ecmwf_data_91(:,2)/24));
    n_levels_91 = round(size(ecmwf_data_91,1)/n_time_91);
  
    if n_time_91*n_levels_91~=size(ecmwf_data_91,1)
        disp('ERROR: Mix-up of ECMWF data...')
        return
    end
    
    % Make a matrix out of the data vectors with 91 levels:
    time91 = reshape(ecmwf_data_91(:,1)+ecmwf_data_91(:,2)/24,n_levels_91,n_time_91);
    p91 = reshape(ecmwf_data_91(:,3),n_levels_91,n_time_91);
    T91 = reshape(ecmwf_data_91(:,4),n_levels_91,n_time_91);
    z91 = reshape(ecmwf_data_91(:,5),n_levels_91,n_time_91);
    vmr91 = reshape(ecmwf_data_91(:,6),n_levels_91,n_time_91);
    
    % Interpolate to constant pressure levels (91)
    p_i = power(10,nanmean(log10(p91),2));
    T_i91 = NaN(size(T91));
    z_i91 = NaN(size(z91));
    vmr_i91 = NaN(size(vmr91));
    for i = 1:size(T91,2)
        T_i91(:,i) = interp1(log10(p91(:,i)),T91(:,i),log10(p_i));
        z_i91(:,i) = interp1(log10(p91(:,i)),z91(:,i),log10(p_i));
        vmr_i91(:,i) = interp1(log10(p91(:,i)),vmr91(:,i),log10(p_i));
    end
    
    % Determine number of levels and times with 137 levels:
    n_time_137 = length(unique(ecmwf_data_137(:,1)+ecmwf_data_137(:,2)/24));
    n_levels_137 = round(size(ecmwf_data_137,1)/n_time_137);
  
    if n_time_137*n_levels_137~=size(ecmwf_data_137,1)
        disp('ERROR: Mix-up of ECMWF data...')
        return
    end
    
    % Make a matrix out of the data vectors with 137 levels:
    time137 = reshape(ecmwf_data_137(:,1)+ecmwf_data_137(:,2)/24,n_levels_137,n_time_137);
    p137 = reshape(ecmwf_data_137(:,3),n_levels_137,n_time_137);
    T137 = reshape(ecmwf_data_137(:,4),n_levels_137,n_time_137);
    z137 = reshape(ecmwf_data_137(:,5),n_levels_137,n_time_137);
    vmr137 = reshape(ecmwf_data_137(:,6),n_levels_137,n_time_137);
        
    % Interpolate to constant pressure levels (91)
    T_i137 = NaN(91,size(T137,2));
    z_i137 = NaN(91,size(z137,2));
    vmr_i137 = NaN(91,size(vmr137,2));
    for i = 1:size(T137,2)
        T_i137(:,i) = interp1(log10(p137(:,i)),T137(:,i),log10(p_i));
        z_i137(:,i) = interp1(log10(p137(:,i)),z137(:,i),log10(p_i));
        vmr_i137(:,i) = interp1(log10(p137(:,i)),vmr137(:,i),log10(p_i));
    end
        

    % Put 91 and 137 together
    time = [time91 time137(1:91,:)];
    T_i = [T_i91 T_i137];
    z_i = [z_i91 z_i137];
    vmr_i = [vmr_i91 vmr_i137];
    
    % Take the average in the given time frame
    time_ind = (time(1,:)>=starttime & time(1,:)<=stoptime);
    p = p_i;
    T = nanmean(T_i(:,time_ind),2);
    z = nanmean(z_i(:,time_ind),2);
    vmr = nanmean(vmr_i(:,time_ind),2);
    ptz = flipud([p T z]);
    vmr = flipud(vmr);
    
    ptz = ptz(~isnan(ptz(:,1)) & ~isnan(ptz(:,2)) & ~isnan(ptz(:,3)),:);
    vmr = vmr(~isnan(ptz(:,1)) & ~isnan(ptz(:,2)) & ~isnan(ptz(:,3)));
    
    if do_cira
        % get CIRA data for heights above 80 km
        [T_cira, D_cira] = get_cira86(mean([starttime stoptime]), lat);
        
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
        
        vmr=[vmr; nan(length(ind),1)];
    end    
    
else
    
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
    vmr = reshape(ecmwf_data(:,6),n_levels,n_time);
    
    % Interpolate to constant pressure levels
    p_i = power(10,nanmean(log10(p),2));
    T_i = NaN(size(T));
    z_i = NaN(size(z));
    vmr_i = NaN(size(z));
    for i = 1:size(T,2)
        T_i(:,i) = interp1(log10(p(:,i)),T(:,i),log10(p_i));
        z_i(:,i) = interp1(log10(p(:,i)),z(:,i),log10(p_i));
        vmr_i(:,i) = interp1(log10(p(:,i)),vmr(:,i),log10(p_i));
    end
    
    % Take the average in the given time frame
    time_ind = (time(1,:)>=starttime & time(1,:)<=stoptime);
    p = p_i;
    T = nanmean(T_i(:,time_ind),2);
    z = nanmean(z_i(:,time_ind),2);
    vmr = nanmean(vmr_i(:,time_ind),2);
    ptz = flipud([p T z]);
    vmr = flipud(vmr);
    
    ptz = ptz(~isnan(ptz(:,1)) & ~isnan(ptz(:,2)) & ~isnan(ptz(:,3)),:);
    vmr = vmr(~isnan(ptz(:,1)) & ~isnan(ptz(:,2)) & ~isnan(ptz(:,3)));
    
    
    if do_cira
        % get CIRA data for heights above 80 km
        [T_cira, D_cira] = get_cira86(mean([starttime stoptime]), lat);
        
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
        
        vmr=[vmr; nan(length(ind),1)];
    end
    
end