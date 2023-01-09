% insert_profile(out,version,mission)
% 
% Speichert Retrieval in mysql Datenbank
% out: Retrieval output
% version: specify the retrieval version number and make comments in the version table
% mission: [1: bern, 2: sodankyla, 3: zimmerwald]
function prof = insert_profile_into_mysql(out,version,mission)

% return, if the mission is not 1,2 or 3
if isempty(find(mission==[1 2 3 4]))==1
    disp('invalid mission!')
    return
end

% if aos_quality_check(out)>0
% if out.char.G.niter>2 | sum((out.spectrum(3,:)-out.spectrum(4,:)).^2)>1
% 	disp('pfui dr tüüüfu!');
% 	return
% end

% Datenbank parameter
db='MIAWARA';
% db='SWARA';

if strcmp(db,'SWARA_test')==1
    disp('inserting retrieval into SWARA_test!')
end

if strcmp(db,'SWARA')==1
    disp('inserting retrieval into SWARA!')
end

    

mysql(['use ' db]);

  
% Retrieval parameter
header.lat = out.level1.latitude;
header.lon = out.level1.longitude;
header.alt = out.level1.altitude;
header.time_min = datestr(out.level1.min_time,31);
header.time_max = datestr(out.level1.max_time,31);
header.time = datestr(out.level1.time,31);
header.integration_time = out.level1.tint;
header.zenith_angle = 0;
header.version = version;
header.mission = mission;
    
    
% Header zusammenstellen & einfügen
head_sql = sql_head_string(header);
[a,b,headerid] = mysql(head_sql);
   

% AVK zusammenstellen und einfügen
[m,n]=size(out.char.avk);
AVK.A = out.char.avk;
AVK.number_of_rows =  m;
AVK.number_of_columns = n;
    
% AVK einfügen
sql_AVK_insert(AVK,headerid);
   
% Profilinfo zusammenstellen
for i = 1:m
      prof.ID = headerid;
      prof.alt = out.profile(i,3);
      prof.p = out.profile(i,1);
      prof.T = out.profile(i,2);
      prof.vmr = out.profile(i,4);
      prof.vmr_smo = out.char.s_smo(i) .* out.profile(i,4);
      prof.vmr_obs = out.char.s_obs(i) .* out.profile(i,4);
      prof.apr = out.profile(i,5);
      prof.aprc = 1 - out.char.measres(i);
      
      % Profilinfo in DB einfügen
      prof_sql = sql_prof_string(prof);
      mysql(prof_sql);
end


function head_sql = sql_head_string(header)
% präpariert den SQL insert string für die Profil Tabelle

% Basis info
head_sql = sprintf('INSERT INTO level2_header SET ');
% Zeit info
head_sql = sprintf('%s time = "%s", ',head_sql, header.time);
head_sql = sprintf('%s time_min = "%s", ',head_sql, header.time_min);
head_sql = sprintf('%s time_max = "%s", ',head_sql, header.time_max);
% integrationszeit
head_sql = sprintf('%s integration_time = "%f", ',head_sql, header.integration_time);
% Geographische Info
head_sql = sprintf('%s latitude = "%.3f", ',head_sql, header.lat);
head_sql = sprintf('%s longitude = "%.3f", ',head_sql, header.lon);
head_sql = sprintf('%s altitude = "%.3f", ',head_sql, header.alt);
% Zenithangle
head_sql = sprintf('%s zenith_angle = "%.3f", ',head_sql, header.zenith_angle);
% mission
head_sql = sprintf('%s mission = "%d", ',head_sql, header.mission);
% Version
head_sql = sprintf('%s version = "%3.1f";', head_sql, header.version);




function prof_sql = sql_prof_string(prof)
% präpariert den SQL insert string für die Profil Tabelle
clear prof_sql
% Basis info
prof_sql = sprintf('INSERT INTO level2_profile SET ');
% headerID
prof_sql = sprintf('%s profile_ID = "%d", ', prof_sql, prof.ID);
% PTZ
prof_sql = sprintf('%s altitude = "%.3f", ', prof_sql, prof.alt);
prof_sql = sprintf('%s pressure = "%.4f", ', prof_sql, prof.p);
prof_sql = sprintf('%s temperature = "%.2f", ', prof_sql, prof.T);
% VMR
prof_sql = sprintf('%s vmr = "%5.9f", ', prof_sql, prof.vmr);
prof_sql = sprintf('%s vmr_error_observation = "%5.9f", ', prof_sql, prof.vmr_obs);
prof_sql = sprintf('%s vmr_error_smooth = "%5.9f", ', prof_sql, prof.vmr_smo);
% apriori
prof_sql = sprintf('%s apriori_vmr = "%5.9f", ', prof_sql, prof.apr);
prof_sql = sprintf('%s apriori_contribution = "%f"; ', prof_sql, prof.aprc);


function  sql_AVK_insert(AVK,header_ID)

column = AVK.number_of_columns;
row = AVK.number_of_rows;
AVK_ID =  header_ID;


for i = 1:row
  AVK_sql_0 = sprintf('INSERT INTO level2_AVK SET profile_ID = "%d", AVK_row = "%d", ',header_ID, i); 

  for k = 1:column
    AVK_sql = sprintf('%s  AVK_col = "%d", AVK = "%.9f";',AVK_sql_0,k,AVK.A(i, k));
    mysql(AVK_sql);
  end
end







