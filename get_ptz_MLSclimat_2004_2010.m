function ptz=get_ptz_MLSclimat_2004_2010(time_num)

% select data
C=mysql('select * from atmospheres.MLS_Aura_T_climatology_2004_2010','mat');

% reshape data
p=reshape(C(:,4),52,365);
T=reshape(C(:,5),52,365);
z=reshape(C(:,6),52,365);

% convert time to day of year
d=floor(time_num)-datenum(year(time_num),1,0);

% assembel ptz
ptz=[p(:,d) T(:,d) z(:,d)];

ptz = sortrows(ptz,-1);

% get CIRA data
[T_cira, D_cira] = get_cira86(time_num, 47);

Z_cira_i=(0:1000:120000)';
T_cira_i=interp1(T_cira(:,1),T_cira(:,2),Z_cira_i);
D_cira_i=exp(interp1(D_cira(:,1),log(D_cira(:,2)),Z_cira_i));

ind=find(isnan(T_cira_i)==0 & isnan(D_cira_i)==0);

Z_cira_i=Z_cira_i(ind);
T_cira_i=T_cira_i(ind);
D_cira_i=D_cira_i(ind);

% find CIRA indices above AuraMLS profile
% ind=find(Z_cira_i>ptz(end,3) & D_cira_i<ptz(end,1));
ind=find(Z_cira_i>ptz(end,3) & Z_cira_i>95000 & D_cira_i<ptz(end,1));

% merge the two profiles
ptz=[ptz; D_cira_i(ind) T_cira_i(ind) Z_cira_i(ind)];
