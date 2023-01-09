
function insert_L2_from_files_into_mysql(filenames,retrieval_version)

% filenames{1} = '/export/data/MIAWARA/retrievals/RETRIEVAL_ARTS2/retrieval_331620.01/retrieval_20070328-1020.mat';
% retrieval_version = 331620.01;

mysql('use MIAWARA');

for i = 1:length(filenames)
    
    disp(['Inserting ' filenames{i}]);
    load(filenames{i});
    
    time = datestr([L2.year L2.month L2.day L2.hour L2.minute L2.second],31);
    time_min = datestr(L2.min_time,31);
    time_max = datestr(L2.max_time,31);
    if (strcmp(L2.M.int_type,'fixed_sigma') || strcmp(L2.M.int_type,'fixed_sigma_low_noise'))
        time_min_supposed = time_min;
        time_max_supposed = time_max;
    else
        time_min_supposed = datestr(datenum(L2.M.starttime,31),31);
        time_max_supposed = datestr(datenum(L2.M.stoptime,31),31);
    end
    integration_time = num2str(L2.tint);
    level1_noise = num2str(L2.tnoise);
    latitude = num2str(46.8770);
    longitude = num2str(7.4650);
    altitude = num2str(907);
    version = num2str(retrieval_version);
    p_lower_ind = find(L2.species1_mr>=0.75,1,'first');
    p_upper_ind = find(L2.species1_mr>=0.75,1,'last');
    if ~isempty(p_lower_ind) && ~isempty(p_upper_ind)
        p_lower_lim = L2.species1_p(p_lower_ind);
        p_upper_lim = L2.species1_p(p_upper_ind);
        if p_lower_lim<p_upper_lim
            tmpvar = p_upper_lim;
            p_upper_lim = p_lower_lim;
            p_lower_lim = tmpvar;
        end
        pressure_lower_limit = num2str(p_lower_lim);
        pressure_upper_limit = num2str(p_upper_lim);
    else
        pressure_lower_limit = '0';
        pressure_upper_limit = '0';
    end
    
    values = ['"' time '","' time_min '","' time_max '","' time_min_supposed '","' time_max_supposed '",' integration_time ',' level1_noise ',' latitude ',' longitude ',' altitude ',' pressure_lower_limit ',' pressure_upper_limit ',' version];
    
    % Insert header data into MySQL
    header_sql = ['INSERT INTO level2_header (time,time_min,time_max,time_min_supposed,time_max_supposed,integration_time,level1_noise,latitude,longitude,altitude,press_lower_limit,press_upper_limit,version,zenith_angle,mission) VALUES(' values ',0,3)'];
    [tmpvar,tmpvar2,prof_ID] = mysql(header_sql);
    
    % Insert AVK's into MySQL
    for j = 1:size(L2.species1_A,1)
        
        AVK_sql_0 = sprintf('INSERT INTO level2_AVK SET profile_ID = "%d", AVK_row = "%d", ',prof_ID, j);
        for k = 1:size(L2.species1_A,2)
            AVK_sql = sprintf('%s  AVK_col = "%d", AVK = "%.9f";',AVK_sql_0,k,full(L2.species1_A(j, k)));
            mysql(AVK_sql);
        end
        
    end
    
    % Insert Profile into MySQL
    for j = 1:size(L2.species1_A,1)
        
        prof_sql = sprintf('INSERT INTO level2_profile SET ');
        prof_sql = sprintf('%s profile_ID = "%d", ', prof_sql, prof_ID);
        prof_sql = sprintf('%s altitude = "%.3f", ', prof_sql, interp1(log10(L2.Q.Z_ATMDATA.GRID1),L2.Q.Z_ATMDATA.DATA,log10(L2.species1_p(j))));
        prof_sql = sprintf('%s pressure = "%.4f", ', prof_sql, L2.species1_p(j));
        prof_sql = sprintf('%s temperature = "%.2f", ', prof_sql, interp1(log10(L2.Q.T.ATMDATA.GRID1),L2.Q.T.ATMDATA.DATA,log10(L2.species1_p(j))));
        prof_sql = sprintf('%s vmr = "%5.9f", ', prof_sql, L2.species1_x(j));
        prof_sql = sprintf('%s vmr_error_observation = "%5.9f", ', prof_sql, L2.species1_eo(j));
        prof_sql = sprintf('%s vmr_error_smooth = "%5.9f", ', prof_sql, L2.species1_es(j));
        prof_sql = sprintf('%s apriori_vmr = "%5.9f", ', prof_sql, L2.species1_xa(j));
        prof_sql = sprintf('%s apriori_contribution = "%f"; ', prof_sql, 1 - L2.species1_mr(j));
        
        mysql(prof_sql);
        
    end
    
end