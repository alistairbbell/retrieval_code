% function plotL2(L2,filename)
%
% This functions plot a retrieval done by ARTS2/Qpack2, based on the
% L2-structure. The optional variable what_to_plot is a string and
% can be either 'normal' (default) or 'all'.
%
function f = saveL2plot(L2,filename)

% filename = '/home/scheiben/test.jpg';

% Bern:
lat = 46.877;
lon = 7.465;

% Try to find a satellite profile
if isfield(L2,'min_time')
    A=get_mls_from_mysql(datestr(L2.min_time,31),datestr(L2.max_time,31),lat,lon,'H2O',3.3);
else
    A=get_mls_from_mysql(datestr(datenum([L2.year L2.month L2.day L2.hour L2.minute L2.second])-1.5,31),datestr(datenum([L2.year L2.month L2.day L2.hour L2.minute L2.second])+1.5,31),lat,lon,'H2O',3.3);
end

f = figure('Visible','Off');
subplot(2,2,1);
plot(L2.species1_x*1e6,L2.species1_p/100);
hold on
plot(L2.species1_x(L2.species1_mr>0.8)*1e6,L2.species1_p(L2.species1_mr>0.8)/100,'LineWidth',1.5);
plot(L2.species1_xa*1e6,L2.species1_p/100,'--');
if ~isempty(A); plot(A(:,2)*1e6,A(:,1)/100,'r'); end
set(gca,'YScale','log','YDir','reverse','FontSize',15);
axis([1 8 0.005 50]);
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
set(gca,'YScale','log','YDir','reverse','FontSize',15);
axis([-0.2 1.2 0.005 50]);
grid on
xlabel('AVK*3 and measurement response [-]');
ylabel('Pressure [hPa]');
title('Averaging kernels and measurement response')

subplot(2,2,3);
if isfield(L2,'Q')
    plot(L2.f/1e9,L2.y-L2.bl,'.','MarkerSize',2);
    hold on
    plot(L2.f/1e9,L2.yf-L2.bl,'r','LineWidth',1);
    set(gca,'XLim',[min(L2.f)/1e9 max(L2.f)/1e9],'FontSize',15);
else
    plot(L2.f/1e9,L2.y-L2.bl,'.');
    hold on
    plot(L2.f/1e9,L2.yf-L2.bl,'r','LineWidth',2);
    set(gca,'XLim',[min(L2.f)/1e9 max(L2.f)/1e9],'FontSize',15);
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
    plot(L2.f/1e9,L2.y-L2.yf);
    set(gca,'XLim',[min(L2.f)/1e9 max(L2.f)/1e9]);
else
    plot(L2.f/1e9,L2.y-L2.yf);
    set(gca,'XLim',[min(L2.f)/1e9 max(L2.f)/1e9]);
end
hold on
plot(L2.f/1e9,smooth(L2.y-L2.yf,250),'r','LineWidth',1.5);
grid on
xlabel('Frequency [GHz]');
ylabel('Residuals [K]');
set(gca,'YLim',[-min(L2.Y.TNOISE)*4 min(L2.Y.TNOISE)*4],'FontSize',15);
title('Residuals between measurements and forward model')

save_figure(gcf,filename(1:end-4),7,13,17.5,'jpg',80)

close(f);
