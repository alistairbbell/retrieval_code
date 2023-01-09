% function plotL2(L2,what_to_plot)
%
% This functions plot a retrieval done by ARTS2/Qpack2, based on the
% L2-structure. The optional variable what_to_plot is a string and
% can be either 'normal' (default) or 'all'.
%
function f = plotL2(L2,what_to_plot,additional_data)

if ~exist('what_to_plot','var'); what_to_plot = 'normal'; end
if ~exist('additional_data','var'); additional_data = true; end

% Onsala:
%lat = 57.414722;
%lon = 12.0225;

% Bern:
lat = 46.877;
lon = 7.465;

if additional_data
    % Try to find a satellite profile
    if isfield(L2,'min_time')
        A=get_mls_from_mysql(datestr(L2.min_time,31),datestr(L2.max_time,31),lat,lon,'H2O',3.3);
    else
        A=get_mls_from_mysql(datestr(datenum([L2.year L2.month L2.day L2.hour L2.minute L2.second])-1.5,31),datestr(datenum([L2.year L2.month L2.day L2.hour L2.minute L2.second])+1.5,31),lat,lon,'H2O',3.3);
    end
    % Try to find the ECMWF profile
    if isfield(L2,'min_time')
        ecmwf=get_ecmwf_profiles(datestr(L2.min_time,31),datestr(L2.max_time,31),lat,lon,'H2O_VMR');
    else
        ecmwf=get_ecmwf_profiles(datestr(datenum([L2.year L2.month L2.day L2.hour L2.minute L2.second])-1.5,31),datestr(datenum([L2.year L2.month L2.day L2.hour L2.minute L2.second])+1.5,31),lat,lon,'H2O_VMR');
    end
end

f = figure;
subplot(2,2,1);
if strcmp(L2.Q.ABS_SPECIES(1).UNIT,'vmr')
    plot(L2.species1_x*1e6,L2.species1_p/100);
    hold on
    plot(L2.species1_x(L2.species1_mr>0.8)*1e6,L2.species1_p(L2.species1_mr>0.8)/100,'LineWidth',2);
    plot(L2.species1_xa*1e6,L2.species1_p/100,'--');
elseif strcmp(L2.Q.ABS_SPECIES(1).UNIT,'rel')
    apr_x = L2.Q.ABS_SPECIES(1).ATMDATA.DATA;
    apr_p = L2.Q.ABS_SPECIES(1).ATMDATA.GRID1;
    apr_x_i = interp1(log(apr_p),apr_x,log(L2.species1_p));
    plot(L2.species1_x.*apr_x_i*1e6,L2.species1_p/100);
    hold on
    plot(L2.species1_x(L2.species1_mr>0.8).*apr_x_i(L2.species1_mr>0.8)*1e6,L2.species1_p(L2.species1_mr>0.8)/100,'LineWidth',2);
    plot(apr_x_i*1e6,L2.species1_p/100,'--');
else
    error('Retrieval units not understood!')
end
if additional_data
    if ~isempty(A); plot(A(:,2)*1e6,A(:,1)/100,'r'); end
    if ~isempty(ecmwf); plot(ecmwf(:,2)*1e6,ecmwf(:,1)/100,'g'); end
end
set(gca,'YScale','log','YDir','reverse');
axis([2 8 0.005 250]);
grid on
xlabel('H_2O mixing ratio [ppm]');
ylabel('Pressure [hPa]');
if isfield(L2,'min_time')
    title([datestr(L2.min_time,'yyyy-mm-dd HH:MM') ' - ' datestr(L2.max_time,'yyyy-mm-dd HH:MM')]);
else
    title(datestr([L2.year L2.month L2.day L2.hour L2.minute L2.second],31));
end

subplot(2,2,2);
plot(L2.species1_A'*3,L2.species1_p/100);
hold on
plot(L2.species1_mr+0.01,L2.species1_p/100,'k');
set(gca,'YScale','log','YDir','reverse');
axis([-0.2 1.2 0.005 250]);
grid on
xlabel('AVK*3 and measurement response [-]');
ylabel('Pressure [hPa]');
title('Averaging kernels and measurement response')

subplot(2,2,3);
if isfield(L2,'Q')
    if isfield(L2.Q.SENSOR_RESPONSE,'LO') && ~isempty(L2.Q.SENSOR_RESPONSE.LO)
        plot((L2.f+L2.Q.SENSOR_RESPONSE.LO)/1e9,L2.y-L2.bl,'.');
        hold on
        plot((L2.f+L2.Q.SENSOR_RESPONSE.LO)/1e9,L2.yf-L2.bl,'r','LineWidth',2);
        set(gca,'XLim',[min((L2.f+L2.Q.SENSOR_RESPONSE.LO))/1e9 max((L2.f+L2.Q.SENSOR_RESPONSE.LO))/1e9]);
    else
        plot(L2.f/1e9,L2.y-L2.bl,'.');
        hold on
        plot(L2.f/1e9,L2.yf-L2.bl,'r','LineWidth',2);
        set(gca,'XLim',[min(L2.f)/1e9 max(L2.f)/1e9]);
    end
else
    plot(L2.f/1e9,L2.y-L2.bl,'.');
    hold on
    plot(L2.f/1e9,L2.yf-L2.bl,'r','LineWidth',2);
    set(gca,'XLim',[min(L2.f)/1e9 max(L2.f)/1e9]);
end
grid on
xlabel('Frequency [GHz]');
ylabel('Brightness temperature [K]');
if nanmedian(L2.y-L2.bl)>2.2 && nanmedian(L2.y-L2.bl)<2.65
    set(gca,'YLim',[2.2 2.65]);
end
title('Measurements and forward model')

subplot(2,2,4);
if isfield(L2,'Q')
    if isfield(L2.Q.SENSOR_RESPONSE,'LO') && ~isempty(L2.Q.SENSOR_RESPONSE.LO)
        plot((L2.f+L2.Q.SENSOR_RESPONSE.LO)/1e9,L2.y-L2.yf);
        hold on
        plot((L2.f+L2.Q.SENSOR_RESPONSE.LO)/1e9,smooth(L2.y-L2.yf,250),'r');        
        set(gca,'XLim',[min((L2.f+L2.Q.SENSOR_RESPONSE.LO))/1e9 max((L2.f+L2.Q.SENSOR_RESPONSE.LO))/1e9]);
    else
        plot(L2.f/1e9,L2.y-L2.yf);
        hold on
        plot(L2.f/1e9,smooth(L2.y-L2.yf,250),'r');
        set(gca,'XLim',[min(L2.f)/1e9 max(L2.f)/1e9]);
    end
else
    plot((L2.f)/1e9,L2.y-L2.yf);
    hold on
    plot((L2.f)/1e9,smooth(L2.y-L2.yf,250),'r');    
    set(gca,'XLim',[min(L2.f)/1e9 max(L2.f)/1e9]);
end
grid on
xlabel('Frequency [GHz]');
ylabel('Residuals [K]');
set(gca,'YLim',[-L2.tnoise*4 L2.tnoise*4]);
title('Residuals between measurements and forward model')

if strcmp(what_to_plot,'all')
    
    figure;
    subplot(2,2,1);
    surf(L2.Q.ABS_SPECIES(1).GRIDS{1}/100,L2.Q.ABS_SPECIES(1).GRIDS{1}/100,sqrt(full(L2.Q.ABS_SPECIES(1).SX))*1e6,'LineStyle','None','FaceColor','flat');
    set(gca,'YScale','log','YDir','reverse','XScale','log','XDir','reverse');
    view(0,90);
    h_c1 = colorbar;
    ylabel(h_c1,'H_2O VMR standard deviation [ppm]')
    xlabel('Pressure [hPa]');
    ylabel('Pressure [hPa]');
    axis([min(L2.Q.ABS_SPECIES(1).GRIDS{1}/100) max(L2.Q.ABS_SPECIES(1).GRIDS{1}/100) min(L2.Q.ABS_SPECIES(1).GRIDS{1}/100) max(L2.Q.ABS_SPECIES(1).GRIDS{1}/100)])
    title('H_2O covariance matrix (absolute values)');
    
    switch L2.R.H2O_apriori_covariance
        case 'v17'
            SX_rel = read_datafile(fullfile(fileparts(L2.Q.ABS_LINES),'v17.sx.apriori.H2O.aa'),'aomatrix');
        case 'v18'
            SX_rel = read_datafile(fullfile(fileparts(L2.Q.ABS_LINES),'v18.sx.apriori.H2O.aa'),'aomatrix');
        case 'v22'
            SX_v22_abs = read_datafile(fullfile(fileparts(L2.Q.ABS_LINES),'v22.sx.apriori.H2O.ppm.aa'),'matrix');
            SX_rel{1} = [3;0];
            SX_rel{2}(:,1) = log10(L2.Q.ABS_SPECIES(1).ATMDATA.GRID1);
            SX_rel{2}(:,2) = interp1(SX_v22_abs(:,1),SX_v22_abs(:,2),L2.Q.ABS_SPECIES(1).ATMDATA.GRID1);
            SX_rel{2}(:,3) = 0.25;
            rel = SX_rel{2}(:,2)./L2.Q.ABS_SPECIES(1).ATMDATA.DATA;
            rel(rel<0.15 | rel>4) = 0.15;
            rel(rel>0.65) = 0.65;
            rel(isnan(rel)) = 0.65;
            rel(1) = rel(2);
            rel=smooth(rel,10);
            SX_rel{2}(:,2) = rel;
            clear rel SX_v22_abs
        otherwise
            SX_rel = [];
    end
    if ~isempty(SX_rel)
        switch SX_rel{1}(1)
            case 0
                corr_fun = 'drc';
            case 1
                corr_fun = 'linn';
            case 2
                corr_fun = 'exp';
            case 3
                corr_fun = 'gau';
            otherwise
                error('Correlation function of the H2O covariance matrix not understood!')
        end
        SX_rel_full = covmat1d_from_cfun( L2.Q.ABS_SPECIES(1).GRIDS{1}, [power(10,SX_rel{2}(:,1)) SX_rel{2}(:,2)], corr_fun, [power(10,SX_rel{2}(:,1)) SX_rel{2}(:,3)], SX_rel{1}(2), @log10);
        clear SX_rel
        subplot(2,2,2);
        surf(L2.Q.ABS_SPECIES(1).GRIDS{1}/100,L2.Q.ABS_SPECIES(1).GRIDS{1}/100,sqrt(full(SX_rel_full))*100,'LineStyle','None','FaceColor','flat');
        set(gca,'YScale','log','YDir','reverse','XScale','log','XDir','reverse');
        view(0,90);
        h_c1 = colorbar;
        ylabel(h_c1,'H_2O VMR standard deviation [%]')
        xlabel('Pressure [hPa]');
        ylabel('Pressure [hPa]');
        axis([min(L2.Q.ABS_SPECIES(1).GRIDS{1}/100) max(L2.Q.ABS_SPECIES(1).GRIDS{1}/100) min(L2.Q.ABS_SPECIES(1).GRIDS{1}/100) max(L2.Q.ABS_SPECIES(1).GRIDS{1}/100)])
        title('H_2O covariance matrix (relative values)');
    end
    
    subplot(2,2,3)
    plot(L2.Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.grids{1}/1e3,L2.Q.SENSOR_RESPONSE.BACKEND_CHANNEL_RESPONSE{1}.data);
    xlabel('Frequency [kHz]')
    ylabel('Normalized channel response [-]')
    title('Acqiris FFT channel response')
    
    if isfield(L2.Q.SENSOR_RESPONSE,'LO') && ~isempty(L2.Q.SENSOR_RESPONSE.LO)
        subplot(2,2,4)
        plot((L2.Q.SENSOR_RESPONSE.SIDEBAND_RESPONSE.grids{1}+L2.Q.SENSOR_RESPONSE.LO)/1e9,L2.Q.SENSOR_RESPONSE.SIDEBAND_RESPONSE.data);
        xlabel('Frequency [GHz]')
        ylabel('Normalized sideband/mixer response [-]')
        title('MIAWARA sideband/mixer response')
    end
    
    figure;
    plot(L2.species1_e./L2.species1_xa*100,L2.species1_p/100,'LineWidth',2);
    hold on
    plot(L2.species1_eo./L2.species1_xa*100,L2.species1_p/100,'--');
    plot(L2.species1_es./L2.species1_xa*100,L2.species1_p/100,'-.');
    set(gca,'YScale','log','YDir','reverse');
    grid on
    axis([0 75 0.01 100])
    legend('Total error','Observational error','Smoothing error','Location','SouthEast')
    xlabel('Error [%]');
    ylabel('Pressure [hPa]');
    title('Error estimates')
    
end
