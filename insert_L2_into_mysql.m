
function insert_L2_into_mysql(L2,retrieval_version)

mysql('use MIAWARA');

time = datestr([L2.year L2.month L2.day L2.hour L2.minute L2.second],31);
time_min = datestr(L2.min_time,31);
time_max = datestr(L2.max_time,31);
if isfield(L2,'min_time_supposed') && isfield(L2,'max_time_supposed') && ~isempty(L2.min_time_supposed) && ~isempty(L2.max_time_supposed)
    time_min_supposed = datestr(L2.min_time_supposed,31);
    time_max_supposed = datestr(L2.max_time_supposed,31);
else
    time_min_supposed = time_min;
    time_max_supposed = time_max;
end
integration_time = num2str(L2.tint);
level1_noise = num2str(min(L2.Y.TNOISE));
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
[tmpvar1,tmpvar2,prof_ID] = mysql(header_sql);

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

