load('/home/straub/PostDoc/TEMPERA/data_processing_v1/retrieval/Inversion_FFT/retrieval_data/mls_climat')
load('/home/straub/PostDoc/TEMPERA/data_processing_v1/retrieval/Inversion_FFT/retrieval_data/radiosonde_climat')
arts_xmldata_path = atmlab( 'ARTS_XMLDATA_PATH' );
B=gf_artsxml( fullfile( arts_xmldata_path, 'climatology','cira86', 'cira86.t.xml' ));

p_apriori=z2p_simple(0:1000:100000)/100;
for k=1:12
    ind_T_mls=find(~isnan(mls_climat.T(:,k)));
    mls_p=flipud(mls_climat.p(ind_T_mls,k)/100);
    mls_T=flipud(mls_climat.T(ind_T_mls,k));
    ind_T_rs=find(~isnan(radiosonde_climat.T(:,k)));
    rs_p=radiosonde_climat.p(ind_T_rs,k);
    rs_T=radiosonde_climat.T(ind_T_rs,k);
    [p_min_rs ind_min_rs]=min(rs_p);
    [p_max_mls ind_max_mls]=max(mls_p);
    
    p_overlap=[p_max_mls p_min_rs]
    ind_p_overlap_mls=find(mls_p<p_overlap(1)&mls_p>p_overlap(2));
    ind_p_overlap_rs=find(rs_p<p_overlap(1)&rs_p>p_overlap(2));
    fact=linspace(0,1,length(ind_p_overlap_mls)-1);
    
    T_radiosonde_interp=interp1(log10(rs_p),...
        rs_T, log10(mls_p(ind_p_overlap_mls)));
    
    p_apriori1=[rs_p(1:ind_p_overlap_rs(1)); mls_p(ind_p_overlap_mls(1:end-1)); mls_p(ind_p_overlap_mls(end):end)];
    
    T_apriori(:,k)=interp1(log10(p_apriori1),...
        [rs_T(1:ind_p_overlap_rs(1)); mls_T(ind_p_overlap_mls(1:end-1)).*fact'+T_radiosonde_interp(1:end-1).*flipud(fact'); mls_T(ind_p_overlap_mls(end):end)],log10(p_apriori));
end

p_apriori=p_apriori(2:end-5)'*100;

ind = find(B.GRID1<p_apriori(end))
load('/home/straub/PostDoc/TEMPERA/data_processing_v1/retrieval/matlab_BAS/WACCM_Troll_CO_climat.mat')

T=[(T_apriori(:,1)+T_apriori(:,12))/2 T_apriori (T_apriori(:,1)+T_apriori(:,12))/2];
for k=1:14,apriori_T(:,1,1,k)=[T(2:end-5,k); B.DATA(ind,27,1,k) ];end
apriori_p=[p_apriori; B.GRID1(ind) ];

G.NAME='Temperature climatologie over Bern';
G.SOURCE='Radiosonde climatology Payern produced by O. Staehli and MLS climatology in mysql database for 45N';
G.DATA_NAME='T';
G.DATA_UNIT='K';
G.DATA=apriori_T;
G.GRID1=apriori_p;
G.GRID2=47;
G.GRID3=7;

save('/home/straub/PostDoc/TEMPERA/data_processing_v1/retrieval/Inversion_FFT/retrieval_data/T_apriori','G')


