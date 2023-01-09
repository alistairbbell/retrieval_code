function msm = integrate_level1_dual_combo(start_date_str,end_date_str,mission)
% takes level1 data from both polarizations and combines them 


%msm=get_level1_acq_from_NetCDF_certain_time_period_all_dir_and_pol(start_date_str,end_date_str);
[msm1, header]=get_level1_acq_from_mysql_certain_time_periode_dual(start_date_str,end_date_str,1 ,mission);
[msm2, header]=get_level1_acq_from_mysql_certain_time_periode_dual(start_date_str,end_date_str,2 ,mission);
%[msm1]=get_level1_acq_from_NetCDF_certain_time_period(start_date_str,end_date_str,1);
%[msm2]=get_level1_acq_from_NetCDF_certain_time_period(start_date_str,end_date_str,2);



if isempty(msm1.time) && isempty(msm2.time)
disp(sprintf('/--------------------------------------------------------------------\\'))
disp(sprintf('|                   both polarizations are missing                    |'))
disp(sprintf('\\--------------------------------------------------------------------/'))
 
msm.time=msm1.time;
msm.max_time=msm1.max_time;
msm.min_time=msm1.min_time;
msm.altitude=msm1.altitude;
msm.latitude=msm1.latitude;
msm.longitude=msm1.longitude;
msm.tint=msm1.tint;
msm.f=msm1.f;
msm.y=[];
msm.sigma=[];
msm.y1=msm1.y;
msm.sigma1=msm1.sigma;
msm.y2=msm2.y;
msm.sigma2=msm2.sigma;
msm.tau1=msm1.tau;
msm.tau2=msm2.tau;
msm.trec1=msm1.trec;
msm.trec2=msm2.trec;
msm.mirror_ele1=msm1.mirror_elevation;
msm.mirror_ele2=msm2.mirror_elevation;
msm.mirror_elevation=nanmean([msm1.mirror_elevation msm2.mirror_elevation]);
else
    if isempty(msm1.time) || isempty(msm2.time)
        disp(sprintf('/--------------------------------------------------------------------\\'))
        disp(sprintf('|                    one polarization (partially) missing             |'))
        disp(sprintf('\\--------------------------------------------------------------------/'))
        if isempty(msm1.time)
            y=msm2.y;
            sigma=msm2.sigma;
            msm.f=msm2.f;
            msm.time=msm2.time;
            msm.max_time=msm2.max_time;
            msm.min_time=msm2.min_time;
            msm.altitude=msm2.altitude;
            msm.latitude=msm2.latitude;
            msm.longitude=msm2.longitude;
            msm.tint=msm2.tint;
        elseif isempty(msm2.time)
            y=msm1.y;
            sigma=msm1.sigma;
            msm.f=msm1.f;
            msm.time=msm1.time;
            msm.max_time=msm1.max_time;
            msm.min_time=msm1.min_time;
            msm.altitude=msm1.altitude;
            msm.latitude=msm1.latitude;
            msm.longitude=msm1.longitude;
            msm.tint=msm1.tint;
        end
    else
        y=(1/msm1.sigma^2 + 1/msm2.sigma^2)^-1*(msm1.y/msm1.sigma^2+msm2.y/msm2.sigma^2);
        sigma=(1/msm1.sigma^2+1/msm2.sigma^2)^-0.5;
            msm.f=msm1.f;
            msm.time=msm1.time;
            msm.max_time=msm1.max_time;
            msm.min_time=msm1.min_time;
            msm.altitude=msm1.altitude;
            msm.latitude=msm1.latitude;
            msm.longitude=msm1.longitude;
            msm.tint=msm1.tint;
    end
%========== find sigma from signal ===============
ind = find(msm.f>=22260020487.090274810791 & msm.f<=22264964636.513458251953); % corresponds to spectrometer channels 10648 and 10810
m = length(ind);
[p,s,mu]=polyfit(1:m,y(ind),1);
level1_sigma=y(ind)-polyval(p,1:m,[],mu);
sigma_sig=std(level1_sigma);



msm.y=y;
msm.sigma_calc=sigma;
msm.sigma=sigma_sig;
msm.y1=msm1.y;
msm.sigma1=msm1.sigma;
msm.y2=msm2.y;
msm.sigma2=msm2.sigma;
msm.tau1=msm1.tau;
msm.tau2=msm2.tau;
msm.trec1=msm1.trec;
msm.trec2=msm2.trec;
msm.a1=msm1.a;
msm.a2=msm2.a;
msm.A1=msm1.A;
msm.A2=msm2.A;
msm.mirror_ele1=msm1.mirror_elevation;
msm.mirror_ele2=msm2.mirror_elevation;
msm.mirror_elevation=nanmean([msm1.mirror_elevation; msm2.mirror_elevation]);
end
