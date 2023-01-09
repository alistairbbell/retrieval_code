% function ptz = get_ptz_ECMWF_temporal_interpolation(start_time_str, stop_time_str, lat, lon)
%
% 2003-03-13 bd
% History:
% 2005-04-18 ah: Wechsel auf ECMWF ptz.
%
function ptz = get_ptz_ECMWF_tub_temporal_interpolation(msm, extra_info, Y)

ptz = [];

%A and B coefficients needed for ptz conversion
temp = csvread( extra_info.ecmwf_coeff,1, 1,[1,1,138,2] ) ;
A = temp(:,1);
B = temp(:,2);

for i = 1:length(msm)
    disp('msm i')
    disp(length(msm))
    %find datestring
    disp('msm time')
    datetime_obs = datetime(1970,1,1) + days(msm.time);
    disp(datetime_obs)

    datetime_ecmwf = datetime(msm(i).min_time(1)*(60*60*24),'ConvertFrom','epochtime');
    disp(datetime_ecmwf)
    datetime_ecmwf_plusone = datetime_ecmwf+days(1);
    
    datestring_ecmwf = strcat(num2str(year(datetime_ecmwf),'%02.f'), num2str(month(datetime_ecmwf),'%02.f'),num2str(day(datetime_ecmwf),'%02.f'));
    datestring_ecmwf_plusone = strcat(num2str(year(datetime_ecmwf_plusone),'%02.f'), num2str(month(datetime_ecmwf_plusone),'%02.f'),num2str(day(datetime_ecmwf_plusone),'%02.f'));
    
    yearstring = num2str(year(datetime_ecmwf),'%02.f');
    %disp(yearstring)
    
    % Find closest grid point
    lat_grid = -90:1.125:90;
    [tmpvar,lat_grid_ind] = min(abs(Y.LATITUDE-lat_grid));
    

    lat_grid_pt = lat_grid(lat_grid_ind);
    lon_grid=-180:1.125:178.875;
    [tmpvar,lon_grid_ind] = min(abs(Y.LONGITUDE-lon_grid));
    lon_grid_pt = lon_grid(lon_grid_ind);
    
    file_ecmwf = join([extra_info.ecmwf_dir, extra_info.ecmwf_filename, datestring_ecmwf , '.nc']);
    disp(file_ecmwf) ;
    
    file_ecmwf_plusone = join([extra_info.ecmwf_dir, extra_info.ecmwf_filename, datestring_ecmwf_plusone , '.nc']);
    disp(file_ecmwf_plusone) ;
    
    if exist(file_ecmwf)
        disp('Found ecmwf file')
        lat              = ncread(file_ecmwf,'/lat');
        lon              = ncread(file_ecmwf,'/lon');
        time             = ncread(file_ecmwf,'/time');
        pressure         = ncread(file_ecmwf,'/pressure');
        T         = ncread(file_ecmwf,'/temperature');
        q         = ncread(file_ecmwf,'/specific_humidity');
        Phi_Surf     = ncread(file_ecmwf,'/geopotential');
    else
        fprintf('ECMWF file %s not found.\n',file_ecmwf);
    end
    
    


    if exist(file_ecmwf_plusone)
        disp('Found ecmwf file')
        temp_time = ncread(file_ecmwf_plusone,'/time');
        temp_pres = ncread(file_ecmwf_plusone,'/pressure');
        temp_T = ncread(file_ecmwf_plusone,'/temperature');
        temp_q = ncread(file_ecmwf_plusone,'/specific_humidity');
        temp_Phi_Surf     = ncread(file_ecmwf_plusone,'/geopotential');
    
        time(end+1)       = temp_time(1);
        Phi_Surf(end+1)   = temp_Phi_Surf(1);
        pressure(end+1, :)  = temp_pres(1,:);
        T(end+1, :)         = temp_T(1,:);
        q(end+1, :)         = temp_q(1,:);
    else
        fprintf('ECMWF file %s not found.\n',file_ecmwf_plusone);
    end
    
    %disp('ecmwf_times')
    disp(length(time))
    arr_len = length(A);
    
    A_h(1:arr_len-1)=(1/2).*(A(2:arr_len) + A(1:arr_len-1));
    B_h(1:arr_len-1)=(1/2).*(B(2:arr_len) + B(1:arr_len-1));
    T_h(1:arr_len-1)=(1/2).*(T(2:arr_len) + T(1:arr_len-1));
    
    
    for i = 1:length(time)
        %fprintf('Surface pressure:  %f\n',  pressure(i,end));
        [z(i,:), p(i,:)] = ecmwf_zp_calc( A, B, T(i,:)', pressure(i,end), Phi_Surf(i), q(i,:)' );
    end
        
    if isempty(time)
        disp('No ECMWF data found!')
        return
    end
    
    % Determine number of levels and times
    n_time = length(unique(time(:)));
    n_levels = round(length(pressure(1,:)));
    
    fprintf('N time:  %f\n', n_time);
    fprintf('N levels: %f\n', n_levels);
    
    % Interpolate to constant pressure levels
    %disp(size(T))
    p_i = repmat(power(10,mean(log10(p),1, 'omitnan')), length(time),1);
    T_i = NaN(size(T));
    z_i = NaN(size(z));
    for i = 1:n_time
        try
            T_i(i,:) = interp1(log10(p(i,:)),T(i,:),log10(p_i(i,:)));
            z_i(i,:) = interp1(log10(p(i,:)),z(i,:),log10(p_i(i,:)));
        catch ME
            if (strcmp(ME.identifier,'MATLAB:griddedInterpolant:NonUniqueCompVecsPtsErrId'))
                msg = ['Pressure does not vary at altitude of ', ...
                num2str(z(1,i)),' m asl '];
                disp(['error_caught',  int2str(i)])
                T_i(:,i) = T(:,i);
                z_i(:,i) = z(:,i);
            end
        end
    end
    
    % Interpolate to 1-hourly time vector and then take the average
    formatin =  'yyyymmdd';
    disp([datestring_ecmwf,' '  , datestring_ecmwf_plusone])
    disp(datenum(datestr(datenum(datestring_ecmwf,formatin),'yyyymmddHH'),'yyyymmddHH')- datenum(1970,1,1)) 
    time_i = datenum(datestr(datenum(datestring_ecmwf,formatin),'yyyymmddHH'),'yyyymmddHH'):1/24:datenum(datestr(datenum(datestring_ecmwf_plusone,formatin),'yyyymmddHH'),'yyyymmddHH');
    
    for i = 1:length(time_i)
        time_i(i) = (time_i(i)-datenum(1970,1,1))*24*60*60;
    end
    
    time = double(time);
    %time = repmat(time,1,length(p_i(1,:)));
    
    time_i_2d = repmat(time_i,length(p_i(1,:)),1)';
    p_i_2d = repmat(p_i(1,:),length(time_i(1,:)),1);
    
    posixtime_obs = convertTo(datetime_obs, 'posixtime');
    
    disp('Try Interpolation')

    T = interp1(time,T_i,posixtime_obs);
    disp('done1')

    z =  interp1(time, z_i, posixtime_obs);    
    disp('done2')
    
    % disp('time and Z')
    % disp(time_i_2d)
    % disp(size(z_i2))
    % z = mean(z_i2,1, 'omitnan');

    inx_by_alt = find(z>extra_info.proxy_surf_height)
    

    disp([class(T) class(z) class(p)])
    ptz_temp = flipud([p_i(1,inx_by_alt)' T(inx_by_alt)' z(inx_by_alt)']);
    disp('done 3')
    %T = interp2(time,log(p_i(1,:)),T_i',time_i_2d, log(p_i_2d));
    %T_i2 = interp2(time,log(p_i(1,:)),T_i',time_i_2d, log(p_i_2d));
    %z_i2 = interp2(time,log(p_i(1,:)),z_i',time_i_2d, log(p_i_2d));
    %T = mean(T_i2,1, 'omitnan');
    %z = mean(z_i2,1, 'omitnan');
    
    %Replace with time interpolation!!
    
    
    %     disp('dims')
    %     disp(size(p))
    %     disp(size(T))
    %     disp(size(z))
    %ptz_temp = flipud(ptz);%flipud([p T z]);
    


    ptz_temp = ptz_temp(~isnan(ptz_temp(:,1)) & ~isnan(ptz_temp(:,2)) & ~isnan(ptz_temp(:,3)),:);
    disp('ptz temp done')

    if isempty(ptz_temp(i))
       return
    end

    % get CIRA data for heights above 80 km
    disp('CIRA func input')
    disp(mean([datenum(datestring_ecmwf,formatin) datenum(datestring_ecmwf_plusone,formatin)]))
    disp(Y.LATITUDE)
    [T_cira, D_cira] = get_cira86(mean([datenum(datestring_ecmwf,formatin) datenum(datestring_ecmwf_plusone,formatin)]), Y.LATITUDE);
    
    disp('T_cira')
    disp(T_cira)
    
    Z_cira_i=(0:1000:120000)';
    T_cira_i=interp1(T_cira(:,1),T_cira(:,2),Z_cira_i);
    D_cira_i=exp(interp1(D_cira(:,1),log(D_cira(:,2)),Z_cira_i));
       
    ind=find(isnan(T_cira_i)==0 & isnan(D_cira_i)==0);
    
    Z_cira_i=Z_cira_i(ind);
    T_cira_i=T_cira_i(ind);
    D_cira_i=D_cira_i(ind);
    
    % find CIRA indices above ECMWF profile
    % ind=find(Z_cira_i>ptz(end,3) & D_cira_i<ptz(end,1));
    ind=find(Z_cira_i>max(ptz_temp(:,3)) & Z_cira_i>95000 & D_cira_i<min(ptz_temp(:,1)));
    
    % merge the two profiles
    ptz_temp= [ptz_temp; D_cira_i(ind) T_cira_i(ind) Z_cira_i(ind)];
    
    % get USSTD data for the lowest heights
    ptz_usstd_raw = [[101325 22632.1 5474.89 868.019 110.906 66.9389 3.95642 0.3734 0.01]',[288.15 216.65 216.65 228.65 270.65 270.65 214.65 273.15-86.2 273.15-86.2]',[0 11000 20000 32000 47000 51000 71000 84852 100000]'];
    clear ptz_usstd
    ptz_usstd(:,1) = power(10,linspace(log10(ptz_usstd_raw(1,1)),log10(ptz_usstd_raw(end,1)),100))';
    ptz_usstd(:,2) = interp1(log10(ptz_usstd_raw(:,1)),ptz_usstd_raw(:,2),log10(ptz_usstd(:,1)));
    ptz_usstd(:,3) = interp1(log10(ptz_usstd_raw(:,1)),ptz_usstd_raw(:,3),log10(ptz_usstd(:,1)));
    ptz_usstd(:,2) = smoothdata(ptz_usstd(:,2),5);
    ptz_usstd(:,3) = smoothdata(ptz_usstd(:,3),5);
    
    % find USSTD indices below ECMWF profile
    ind=find(ptz_usstd(:,3)<min(ptz_temp(:,3)) & ptz_usstd(:,1)>max(ptz_temp(:,1)));

    % merge the two profiles
    ptz_temp = [ptz_usstd(ind,:); ptz_temp];
    
    format long
    ptz = ptz_temp;
end
