% insert_tP_level2_into_mysql(L2)
% 
% Transfers retrieval output into mysql-database
% L2: Qpack2-Retrieval output
%
% mysql-database: MIAWARA
% mysql-tables:   - level2_header_tp_new     1D-header-data
%                 - level2_profile_tP_new    2D-profile-data
%                 - level2_AVK_tP_new        AVK
%                 - level2_char              characterisation
%                 - level2_avk_params_new    2D-avk parameters (FWHM,...)
%
% Rene Bleisch 2011, rene.bleisch@iap.unibe.ch
function insert_L2_tropo_into_mysql(L2)
% Select DB
mysql('use MIAWARA;');
  
% header
header.lat       = L2.level1.latitude;
header.lon       = L2.level1.longitude;
header.time_min  = datestr(L2.min_time,31);
header.time_max  = datestr(L2.max_time,31);
header.time_mins = datestr(L2.min_time_supposed,31);
header.time_maxs = datestr(L2.max_time_supposed,31);
header.time      = datestr(L2.level1.time,31);
header.timestamp = datestr(now,31);
header.integration_time = L2.M.tint;
header.zenith_angle = 30;
header.sigma     = L2.level1.sigma;
header.tau       = L2.level1.tau;
header.cloud     = L2.R.cloud;
header.prcp      = L2.R.prcp;
header.version    = L2.R.version;

headerid = sql_header_insert(header);

% profile
nz = length(L2.species1_p);
for i = 1:nz
      prof.alt = L2.species1_z(i);
      prof.p   = L2.species1_p(i);
      prof.T   = L2.t_field(i);
      prof.vmr = L2.species1_x(i);
      prof.vmr_smo = L2.char.esrel(i)*L2.species1_x(i);
      prof.vmr_obs = L2.char.eorel(i)*L2.species1_x(i);
      prof.apr     = L2.species1_xa(i);
      prof.aprc    = 1 - L2.species1_mrrel(i);
      
      sql_profile_insert(prof,headerid);
end

% AVK
[m,n]=size(L2.species1_Arel);
AVK.A = L2.species1_Arel;
AVK.number_of_rows =  m;
AVK.number_of_columns = n;

sql_avk_insert(AVK,headerid);

% Characterisation
sql_char_insert(L2.char,headerid);

end




function headerid = sql_header_insert(header)
% inserts header data

head_sql = sprintf('INSERT INTO level2_trop_header SET ');
head_sql = sprintf('%s timestamp = "%s", ', head_sql, header.timestamp);
head_sql = sprintf('%s time = "%s", ',head_sql, header.time);
head_sql = sprintf('%s time_min = "%s", ',head_sql, header.time_min);
head_sql = sprintf('%s time_max = "%s", ',head_sql, header.time_max);
head_sql = sprintf('%s time_min_supposed = "%s", ',head_sql, header.time_mins);
head_sql = sprintf('%s time_max_supposed = "%s", ',head_sql, header.time_maxs);
head_sql = sprintf('%s integration_time = "%f", ',head_sql, header.integration_time);
head_sql = sprintf('%s latitude = "%.3f", ',head_sql, header.lat);
head_sql = sprintf('%s longitude = "%.3f", ',head_sql, header.lon);
%head_sql = sprintf('%s altitude = "%.3f", ',head_sql, header.alt); set to
% 907m by default
head_sql = sprintf('%s zenith_angle = "%.3f", ',head_sql, header.zenith_angle);
head_sql = sprintf('%s T_noise = "%f", ',head_sql, header.sigma);
head_sql = sprintf('%s tau = "%5.9f", ',head_sql, header.tau);
head_sql = sprintf('%s cloud = "%d", ',head_sql, header.cloud);
head_sql = sprintf('%s prcp = "%.3f", ',head_sql, header.prcp);
head_sql = sprintf('%s version = "%3.1f";', head_sql, header.version);

[a,b,headerid] = mysql(head_sql);

end


function prof_sql = sql_profile_insert(prof,headerid)
% inserts profile data

prof_sql = sprintf('INSERT INTO level2_trop_profile SET ');
prof_sql = sprintf('%s profile_ID = "%d", ', prof_sql, headerid);
prof_sql = sprintf('%s altitude = "%.3f", ', prof_sql, prof.alt);
prof_sql = sprintf('%s pressure = "%.2f", ', prof_sql, prof.p);
prof_sql = sprintf('%s temperature = "%.2f", ', prof_sql, prof.T);
prof_sql = sprintf('%s vmr = "%5.9f", ', prof_sql, prof.vmr);
prof_sql = sprintf('%s vmr_error_observation = "%5.9f", ', prof_sql, prof.vmr_obs);
prof_sql = sprintf('%s vmr_error_smooth = "%5.9f", ', prof_sql, prof.vmr_smo);
prof_sql = sprintf('%s apriori_vmr = "%5.9f", ', prof_sql, prof.apr);
prof_sql = sprintf('%s apriori_contribution = "%f"; ', prof_sql, prof.aprc);
mysql(prof_sql);

end

function  sql_avk_insert(AVK,headerid)
% inserts AVK

column = AVK.number_of_columns;
row = AVK.number_of_rows;

for i = 1:row
  AVK_sql_0 = sprintf('INSERT INTO level2_trop_AVK SET profile_ID = "%d", AVK_row = "%d", ',headerid, i); 

  for k = 1:column
    AVK_sql = sprintf('%s  AVK_col = "%d", AVK = "%.9f";',AVK_sql_0,k,AVK.A(i, k));
    mysql(AVK_sql);
  end
end
end

function  sql_char_insert(char,headerid)
% inserts characterisation

% 1D stuff
head_sql = sprintf('INSERT INTO level2_trop_char_1D SET ');
head_sql = sprintf('%s profile_ID= "%d", ', head_sql, headerid);
head_sql = sprintf('%s offset = "%.3f", ', head_sql, char.offset);
head_sql = sprintf('%s slope  = "%.3f", ', head_sql, char.slope);
head_sql = sprintf('%s zMRmax = "%.3f", ', head_sql, char.ap.zMRmax);
head_sql = sprintf('%s MRmax = "%.3f", ', head_sql, char.ap.MRmax);
head_sql = sprintf('%s zMR08b = "%.3f", ', head_sql, char.ap.MR08b);
head_sql = sprintf('%s zMR08t = "%.3f", ', head_sql, char.ap.MR08t);
%head_sql = sprintf('%s rank = "%d", ', head_sql, char.rank);
head_sql = sprintf('%s ds = "%.3f", ', head_sql, char.ds);
head_sql = sprintf('%s H = "%.3f"; ', head_sql, char.H);
mysql(head_sql);

clear head_sql

% 2D stuff
for i=1:length(char.ap.z_nom)
    head_sql = sprintf('INSERT into level2_trop_char_2D SET ');
    head_sql = sprintf('%s profile_ID = "%d", ', head_sql, headerid);
    head_sql = sprintf('%s z_nom = "%.3f", ', head_sql, char.ap.z_nom(i));
    head_sql = sprintf('%s z_avkmax = "%.3f", ', head_sql, char.ap.z_avkmax(i));
    head_sql = sprintf('%s avkmax = "%.9f", ', head_sql, char.ap.avkmax(i));
    head_sql = sprintf('%s z_avkmaxf = "%.3f", ', head_sql, char.ap.z_avkmaxf(i));
    head_sql = sprintf('%s avkmaxf = "%.9f", ', head_sql, char.ap.avkmaxf(i));
    head_sql = sprintf('%s FWHM = "%.3f";', head_sql, char.ap.fwhm(i));
    mysql(head_sql);
    
    clear head_sql
end

end