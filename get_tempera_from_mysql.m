function A=get_tempera_from_mysql(time1,time2,version,take_avk)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if AVK should be taken or not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
    take_avk=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define mysql-database and tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mysql('use TEMPERA;')
header_table='level2_header';
profile_table='level2_profile';
AVK_table='level2_AVK';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get header data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_str1='profile_ID,ROUND(TO_DAYS(time)+1+TIME_TO_SEC(time)/86400,6) AS time,';
var_str2='ROUND(TO_DAYS(time_min)+1+TIME_TO_SEC(time_min)/86400,6) AS time_min,';
var_str3='ROUND(TO_DAYS(time_max)+1+TIME_TO_SEC(time_max)/86400,6) AS time_max,';
var_str4='integration_time,sigma';

cmd=sprintf('select %s%s%s%s from %s where time between "%s" and "%s" and version=%i order by time;'...
    ,var_str1,var_str2,var_str3,var_str4,header_table,time1,time2,version);
header = mysql(cmd,'mat');

if isempty(header)
    disp('no profile in database')
    A=[];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write header data to structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A.id=header(:,1);
A.time=header(:,2);
A.time_min=header(:,3);
A.time_max=header(:,4);
A.integration_time=header(:,5);
A.sigma=header(:,6);    
for k=1:size(header,1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get header data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    var_str1='pressure,altitude,T,T_error_observation,apriori_T';
    var_str2=',apriori_contribution,T_error_smooth';
    
    cmd=sprintf('select %s%s from %s where profile_ID=%i order by pressure;',...
        var_str1,var_str2,profile_table,A.id(k));
    x=mysql(cmd,'mat');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write profile data to structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    A.p(:,k)=x(:,1);
    A.z(:,k)=x(:,2);
    A.T(:,k)=x(:,3);
    A.T_observation_error(:,k)=x(:,4);
    A.T_apriori(:,k)=x(:,5);
    A.T_mr(:,k)=1-x(:,6);
    
    if take_avk==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get AVK data if take_avk
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cmd=sprintf('select AVK_row, AVK_col, AVK from level2_AVK where profile_ID=%i',A.id(k));
        avk_row=mysql(cmd,'mat');
        avk_sparse=sparse(avk_row(:,1),avk_row(:,2),avk_row(:,3));
        X=full(avk_sparse);
        A.avk(k).avk=X;
        
    elseif take_avk==0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % AVK=empty else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A.avk=[];
    end
end
