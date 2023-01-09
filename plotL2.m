% function plotL2(L2)
%
% This functions plot a retrieval done by ARTS2/Qpack2, based on the
% L2-structure.
%
function f = plotL2(L2,font_size_texts)

if ~exist('font_size_texts','var'); font_size_texts = 10; end

%Try to find a satellite profile
try
    A=get_mls_from_mysql(datestr(L2.min_time,31),datestr(L2.max_time,31),L2.Y.LATITUDE,L2.Y.LONGITUDE,'H2O',3.3);
catch
    A = [];
    disp('Could not find MLS H2O data.')
end

% Try to find the ECMWF profile
try
    ecmwf=get_ecmwf_profiles(datestr(L2.min_time,31),datestr(L2.max_time,31),L2.Y.LATITUDE,L2.Y.LONGITUDE,'H2O_vmr');
catch
    ecmwf = [];
    disp('Could not find ECMWF H2O data.')
end

f = figure;
set(f,'Position',[77 151 1040 778]);

p_lim = [0.005 100];
p_ind = find(L2.species1_p/100>=min(p_lim) & L2.species1_p/100<=max(p_lim));
if p_ind(1)==1 && p_ind(end)<length(L2.species1_p)
    p_ind = [p_ind; p_ind(end)+1];
elseif p_ind(1)>1 && p_ind(end)<length(L2.species1_p)
    p_ind = [p_ind(1)-1; p_ind; p_ind(end)+1];
elseif p_ind(1)>1 && p_ind(end)==length(L2.species1_p)
    p_ind = [p_ind(1)-1; p_ind];
end

axes('Position',[0.1*2/3 0.435 0.45*2/3 0.51])
h = NaN(4,1);
h(1) = plot(L2.species1_x*1e6,L2.species1_p/100);
hold on
plot(L2.species1_x(L2.species1_mr>0.8)*1e6,L2.species1_p(L2.species1_mr>0.8)/100,'LineWidth',2);
h(2) = plot(L2.species1_xa*1e6,L2.species1_p/100,'--');
if ~isempty(A); h(3) = plot(A(:,2)*1e6,A(:,1)/100,'r'); end
if ~isempty(ecmwf); h(4) = plot(ecmwf.PV*1e6,ecmwf.pressure/100,'g'); end
set(gca,'YScale','log','YDir','reverse');
axis([1 8 p_lim]);
grid on
xlabel('H_2O mixing ratio [ppm]');
ylabel('Pressure [hPa]');
set(gca,'YTickLabel',{'0.01','0.1','1','10','100'})
title([datestr(L2.min_time,'yyyy-mm-dd HH:MM') ' - ' datestr(L2.max_time,'yyyy-mm-dd HH:MM')]);
if isempty(A) && isempty(ecmwf)
    legend(h(1:2),'Retrieval','A priori','Location','West');
elseif isempty(A) && ~isempty(ecmwf)
    legend(h([1:2 4]),'Retrieval','A priori','ECMWF','Location','West');
elseif ~isempty(A) && isempty(ecmwf)
    legend(h([1:3]),'Retrieval','A priori','Aura MLS','Location','West');
elseif ~isempty(A) && ~isempty(ecmwf)
    legend(h,'Retrieval','A priori','Aura MLS','ECMWF','Location','West');
end

axes('Position',[0.65*2/3 0.585 0.3*2/3 0.36]);
plot(L2.f/1e9,L2.y-L2.bl,'.');
hold on
plot(L2.f/1e9,L2.yf-L2.bl,'r','LineWidth',1.25);
set(gca,'XLim',[min(L2.f)/1e9 max(L2.f)/1e9],'XTickLabel',{''});
grid on
%xlabel('Frequency [GHz]');
ylabel('Tb [K]');
if nanmedian(L2.y-L2.bl)>2.2 && nanmedian(L2.y-L2.bl)<2.65
    set(gca,'YLim',[2.2 2.65]);
end
title('Y, F(Y) and Y-F(Y)');
set(gca,'YLim',[(floor(min(L2.yf-L2.bl)*20)-1)/20 (ceil(max(L2.yf-L2.bl)*20)+1)/20])

axes('Position',[0.65*2/3 0.435 0.3*2/3 0.145]);
plot(L2.f/1e9,L2.y-L2.yf);
hold on
plot(L2.f/1e9,smooth(L2.y-L2.yf,round(length(L2.f)/10)),'g','LineWidth',1.25);
set(gca,'XLim',[min(L2.f)/1e9 max(L2.f)/1e9]);
grid on
xlabel('Frequency [GHz]');
ylabel('Res. [K]');
set(gca,'YLim',std(L2.y-L2.yf)*5*[-1 1]);
y_tick_label = [];
tmpvar = get(gca,'YTick');
for i = 1:length(tmpvar);
    %if i == length(tmpvar)
    %    y_tick_label{i} = '';
    %else
    y_tick_label{i} = num2str(tmpvar(i));
    %end
end
set(gca,'YTickLabel',y_tick_label);

axes('Position',[(0.1)*2/3 0.1 0.82/4*2/3 0.25])
plot(L2.species1_A'*3,L2.species1_p/100,'Color',[0.6 0.6 0.6]);
hold on
plot(L2.species1_mr+0.01,L2.species1_p/100,'k');
[~,tmp_ind] = min(abs(log(L2.species1_p/100)-log(0.1)));
plot(L2.species1_A(tmp_ind,:)'*3,L2.species1_p/100,'b');
text(0.03+max(L2.species1_A(tmp_ind,:)'*3),L2.species1_p(tmp_ind)/100,[num2str(L2.species1_p(tmp_ind)/100,2) ' hPa'],'Color','b','FontSize',font_size_texts)
[~,tmp_ind] = min(abs(log(L2.species1_p/100)-log(1)));
plot(L2.species1_A(tmp_ind,:)'*3,L2.species1_p/100,'b');
text(0.03+max(L2.species1_A(tmp_ind,:)'*3),L2.species1_p(tmp_ind)/100,[num2str(L2.species1_p(tmp_ind)/100,2) ' hPa'],'Color','b','FontSize',font_size_texts)
[~,tmp_ind] = min(abs(log(L2.species1_p/100)-log(10)));
plot(L2.species1_A(tmp_ind,:)'*3,L2.species1_p/100,'b');
text(0.03+max(L2.species1_A(tmp_ind,:)'*3),L2.species1_p(tmp_ind)/100,[num2str(L2.species1_p(tmp_ind)/100,2) ' hPa'],'Color','b','FontSize',font_size_texts)
set(gca,'YScale','log','YDir','reverse');
set(gca,'YTick',[0.01 0.1 1 10 100])
set(gca,'YTickLabel',{'0.01','0.1','1','10','100'})
axis([-0.2 1.3 p_lim]);
grid on
set(gca,'YMinorGrid','off')
xlabel('AVK*3 and MR [-]');
ylabel('Pressure [hPa]');

axes('Position',[(0.1+(0.82/4+0.01)*1)*2/3 0.1 0.82/4*2/3 0.25])
try
    tmpfield = interp1(log(L2.p_grid),L2.z_field,log(L2.species1_p));
    
    tmpvar = fwhm(tmpfield/1000,full(L2.species1_A)); %([50],[50x50]) matrix sparse matrix --> convert to full matrix
   
    plot(tmpvar,L2.species1_p/100);
    set(gca,'YScale','log','YDir','reverse');
    set(gca,'YTick',[0.01 0.1 1 10 100])
    set(gca,'YLim',p_lim);
    set(gca,'XLim',[floor(min(tmpvar(p_ind)))-1 ceil(max(tmpvar(p_ind)))+1]);
    set(gca,'YTickLabel','');
    grid on
    set(gca,'YMinorGrid','off')
    xlabel('Vert. res. [km]');
end

axes('Position',[(0.1+(0.82/4+0.01)*2)*2/3 0.1 0.82/4*2/3 0.25])
AVK_max_p_level = NaN(size(L2.species1_A,1),1);
for i = 1:size(L2.species1_A,1)
    [~,max_ind] = max(L2.species1_A(i,:));
    try
        AVK_max_p_level(i) = L2.species1_p(max_ind);
    end
end
plot(L2.species1_p/100,L2.species1_p/100,'r-');
hold on
plot(AVK_max_p_level/100,L2.species1_p/100,'.-');
set(gca,'YScale','log','YDir','reverse','XScale','log','XDir','reverse');
set(gca,'YTick',[0.01 0.1 1 10 100])
set(gca,'XTick',[0.01 0.1 1 10 100])
set(gca,'YLim',p_lim);
set(gca,'XLim',p_lim);
set(gca,'YTickLabel','');
set(gca,'XTickLabel',{'0.01','0.1','1','10','100'})
grid on
set(gca,'YMinorGrid','off')
set(gca,'XMinorGrid','off')
xlabel('Level of AVK max. [hPa]');

axes('Position',[(0.1+(0.82/4+0.01)*3)*2/3 0.1 0.82/4*2/3 0.25])
factor = 1./L2.species1_xa*100; %%% Relative
% factor = 1e6;    %%% Absolute
plot(L2.species1_eo.*factor,L2.species1_p/100,'r');
hold on
plot(L2.species1_es.*factor,L2.species1_p/100,'g');
plot(L2.species1_e.*factor,L2.species1_p/100,'b');
set(gca,'YScale','log','YDir','reverse');
set(gca,'YTick',[0.01 0.1 1 10 100])
set(gca,'YLim',p_lim);
set(gca,'XLim',[0 (ceil(max([L2.species1_e(p_ind).*factor(p_ind); L2.species1_eo(p_ind).*factor(p_ind); L2.species1_es(p_ind).*factor(p_ind)]*10))+1)/10]);
% set(gca,'XLim',[0 (ceil(max([L2.species1_e(p_ind); L2.species1_eo(p_ind); L2.species1_es(p_ind)]*1e6*10))+1)/10]);
text(2,0.01,'observational','Color','r','FontSize',font_size_texts)
text(2,0.030103,'smoothing','Color','g','FontSize',font_size_texts)
text(2,0.1,'total','Color','b','FontSize',font_size_texts)
set(gca,'YTickLabel','');
grid on
set(gca,'YMinorGrid','off')
%xlabel('Errors [ppm]');
xlabel('Errors [%]');

desc_str = sprintf('Start time: %s',datestr(L2.min_time,31));
desc_str = sprintf('%s\nStop time: %s',desc_str,datestr(L2.max_time,31));
try
    desc_str = sprintf('%s\nActual integration time: %s h',desc_str,num2str(L2.tint/3600));
end
desc_str = sprintf('%s\nNumber of single spectra: %s',desc_str,num2str(size(L2.msm.tau,1)));
try
    desc_str = sprintf('%s\nSpectrometer: %s',desc_str,L2.main_settings.spectrometer);
end
desc_str = sprintf('%s\nBandwidth: %s MHz',desc_str,num2str(L2.main_settings.bandwidth_MHz));
try
    if L2.main_settings.low_noise_int
        desc_str = sprintf('%s\nLow noise integration: Yes',desc_str);
    else
        desc_str = sprintf('%s\nLow noise integration: No',desc_str);
    end
end
try
    if isempty(L2.main_settings.level1_hour_interval)
        desc_str = sprintf('%s\nLevel1 hour interval: Undefined',desc_str);
    else
        desc_str = sprintf('%s\nLevel1 hour interval: [%s]',desc_str,num2str(L2.main_settings.level1_hour_interval));
    end
end
try
    if isempty(L2.main_settings.level1_month)
        desc_str = sprintf('%s\nLevel1 month: Undefined',desc_str);
    else
        desc_str = sprintf('%s\nLevel1 month: [%s]',desc_str,num2str(L2.main_settings.level1_month));
    end
end
try
    if isempty(L2.main_settings.level1_ref_angle)
        desc_str = sprintf('%s\nLevel1 reference angle: Undefined',desc_str);
    else
        desc_str = sprintf('%s\nLevel1 reference angle: [%s] degrees',desc_str,num2str(L2.main_settings.level1_ref_angle));
    end
end
desc_str = sprintf('%s\nLevel1 noise at center: %s K',desc_str,num2str(min(L2.Y.TNOISE)));
desc_str = sprintf('%s\nLevel1 noise at wings: %s K',desc_str,num2str(L2.main_settings.y_noise_at_wings));
desc_str = sprintf('%s\nLevel1 noise, width of gauss shape: %s MHz',desc_str,num2str(L2.main_settings.y_sigma_MHz));
desc_str = sprintf('%s\nPolyfit order: %s',desc_str,num2str(L2.main_settings.polyfit_order));
try
    if L2.main_settings.poly_incl_2nd_iter
        desc_str = sprintf('%s\nPolyfit subtracted before 2nd iteration: Yes',desc_str);
    else
        desc_str = sprintf('%s\nPolyfit subtracted before 2nd iteration: No',desc_str);
    end
end
try
    desc_str = sprintf('%s\nNumber of sinefit coeffs: %s',desc_str,num2str(L2.main_settings.nof_sinefit_coeff));
end
try
    desc_str = sprintf('%s\nSinefit method: %s',desc_str,L2.main_settings.sinefit_method);
end
if ischar(L2.main_settings.ptz_source)
    desc_str = sprintf('%s\nPTZ source: %s',desc_str,strrep(L2.main_settings.ptz_source,'_','\_'));
else
    desc_str = sprintf('%s\nPTZ source: %s',desc_str,'custom');
end
if ischar(L2.main_settings.h2o_apriori)
    desc_str = sprintf('%s\nH2O a priori source: %s',desc_str,strrep(L2.main_settings.h2o_apriori,'_','\_'));
else
    desc_str = sprintf('%s\nH2O a priori source: %s',desc_str,'custom');
end
try
    desc_str = sprintf('%s\nH2O covariance matrix version: %s',desc_str,L2.main_settings.h2o_covmat_vers);
end
try
    desc_str = sprintf('%s\nH2O model: %s',desc_str,strrep(L2.main_settings.h2o_model,'_','\_'));
end
try
    desc_str = sprintf('%s\nH2O line file: %s',desc_str,strrep(L2.main_settings.line_filename,'_','\_'));
end
try
    desc_str = sprintf('%s\nH2O line file format: %s',desc_str,strrep(L2.main_settings.abs_lines_format,'_','\_'));
end
try
    desc_str = sprintf('%s\nH2O line shape: %s',desc_str,strrep(L2.main_settings.line_shape,'_','\_'));
end
try
    if L2.main_settings.line_shape_cutoff==-1
        desc_str = sprintf('%s\nH2O line cut off: No cut off',desc_str);
    else
        desc_str = sprintf('%s\nH2O line cut off: %s Hz',desc_str,num2str(L2.main_settings.line_shape_cutoff));
    end
end
desc_str = sprintf('%s\nH2O line shape factor: %s',desc_str,strrep(L2.main_settings.line_shape_factor,'_','\_'));
try
    desc_str = sprintf('%s\nArtificial instrument altitude: %s m',desc_str,num2str(L2.main_settings.z_platform));
end
try
    if ischar(L2.main_settings.freq_offset)
        desc_str = sprintf('%s\nFrequency offset: %s',desc_str,L2.main_settings.freq_offset);
        desc_str = sprintf('%s\nFrequency offset value: %s kHz',desc_str,num2str(L2.msm.freq_offset_from_pre_retrieval/1000));
    else
        desc_str = sprintf('%s\nFrequency offset: %s kHz',desc_str,num2str(L2.main_settings.freq_offset/1000));
    end
end
desc_str = sprintf('%s\nOEM profile cost: %s',desc_str,num2str(L2.cost_x));
desc_str = sprintf('%s\nOEM measurement cost: %s',desc_str,num2str(L2.cost_y));
desc_str = sprintf('%s\nOEM total cost: %s',desc_str,num2str(L2.cost));
desc_str = sprintf('%s\nOpacity, minimum: %s',desc_str,num2str(min(L2.msm.tau(:,2))));
desc_str = sprintf('%s\nOpacity, average: %s',desc_str,num2str(mean(L2.msm.tau(:,2))));
desc_str = sprintf('%s\nOpacity, maximum: %s',desc_str,num2str(max(L2.msm.tau(:,2))));
desc_str = sprintf('%s\nReceiver temperature, minimum: %s K',desc_str,num2str(min(L2.msm.trec(:,2))));
desc_str = sprintf('%s\nReceiver temperature, average: %s K',desc_str,num2str(mean(L2.msm.trec(:,2))));
desc_str = sprintf('%s\nReceiver temperature, maximum: %s K',desc_str,num2str(max(L2.msm.trec(:,2))));
try
    desc_str = sprintf('%s\nMeasured reference elevations: [%s] degrees',desc_str,num2str(unique(L2.msm.ref_elevation(:,2)')));
end
desc_str = sprintf('%s\nTropospheric correction factor, average: %s',desc_str,num2str(mean(L2.msm.a(:,2))));

h_anno = annotation('textbox',[0.65 0.1 0.33 0.845]);
set(h_anno,'String',desc_str,'FontSize',font_size_texts);
drawnow;
%save_figure(f,'/export/home/miawara/L2_201505',12,20,20,'jpg');

