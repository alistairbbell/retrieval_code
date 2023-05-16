function retrievalPlotting(L2,Y, settings_metadata, savefile)
%%%%% functions to plot and save the retrieval data %%%%
try
    mkdir(savefile)
catch
    'Directory already exists'
end

%save a text file containing important plot information
fn = fieldnames(settings_metadata);
comment = '';
for k=1:numel(fn)
    temp = fn{k};
    comment = append(comment, temp,': ', settings_metadata.(fn{k}), sprintf('\n'));
end

fid = fopen(append(savefile,'plot_info.txt'),'wt');
fprintf(fid, comment);
fclose(fid);


%define the spectra variables to plot
freq = L2.Y.F; %frequency vector
y_sim = L2.yf; % - L2.bl; %retrieval with baseline fit
y_obs = L2.y; % - L2.bl; %original observation
bl = L2.bl; %baseline fit of spectra
res = y_obs - y_sim; %residuals
res_movmean = movmean(res,25);%residual with moving mean
freq_movmean = movmean(freq,25);
Y_meas_err = Y.TNOISE;
A = L2.species1_A; %averaging kernel matrix


%define retrieval quantities
P = L2.species1_p;
z = L2.species1_z;
Xb = L2.species1_xa; % wv background 
X_ret = L2.species1_x;% wv retrieval
mr = L2.species1_mr; %measurement response
T = L2.species1_T; %temperature
Z = L2.species1_z; %altitdude

%create date subfolder (might change later to plot type)
xlims_f = [-.04;.04];

%define struct for freq plots 
PB_freq.xlab = 'Frequency (GHz)';
PB_freq.ylab =  'Brightness Temperature (K)';
PB_freq.logy = false;
PB_freq.ydec = false;
PB_freq.xlims = xlims_f;

%define PB for spectra plots
PB = PB_freq;
PB_freq.logy = false;
PB_freq.ydec = false;
PB.savestring =  append(savefile,'Tb_y_yret');
PB.labels = ["Observation" "Retrieval"];
PB.multi_y = true;
%make plot
plot_basic(freq/1e9-22.235, [(y_obs), (y_sim)], PB);
%clear PB
clear PB

%y error plots 
PB = PB_freq;
PB.savestring =  append(savefile,'Y_noise');
PB.labels = ["Measurement Noise"];
PB.multi_y = false;
PB.xlims = [-.05, .05];
PB.ylims =[0,0.1];
PB.logy = false;
%plot
plot_basic(freq_movmean/1e9-22.235, [Y_meas_err], PB);

PB.xlims = [-.5, .5];
PB.savestring =  append(savefile,'Y_noise_wide');
PB.ylims =[0.0001,100];
PB.logy = true;
plot_basic(freq_movmean/1e9-22.235, [Y_meas_err], PB);


%residual plots 
PB = PB_freq;
PB.savestring =  append(savefile,'Tb_residual_narrow');
PB.labels = ["Residual"];
PB.multi_y = false;
PB.xlims = xlims_f;
PB.ylims = [-.1, 0.1];
%plot
plot_basic(freq_movmean/1e9-22.235, [res_movmean], PB);

%residual plots 
PB = PB_freq;
PB.savestring =  append(savefile,'Tb_residual_wide');
PB.xlims =  [-.5, .5];
PB.labels = ["Residual"];
%plot
plot_basic(freq_movmean/1e9-22.235, [res_movmean], PB);

%Baseline plots
PB = PB_freq;
PB.savestring =  append(savefile,'baseline_narrow');
PB.labels =  ["Baseline"];
%plot
plot_basic(freq/1e9-22.235, [bl ], PB);

%baseline with wide view
PB.savestring =  append(savefile,'baseline_wide');
PB.xlims =  [-.5, .5];
plot_basic(freq/1e9-22.235, [bl ], PB);

%%%   define struct for atmosphere plots %%%
PB_atm.ylab =  'Pressure (Pa)';
PB_atm.logy = true;
PB_atm.ydec = true;
PB_atm.xlims = false;

%make plot for retrieval
clear PB
PB = PB_atm;
PB.xlab = 'H_{2}O (PPMV)';
PB.ylab = 'P (Pa)';
PB.savestring = append(savefile,'wv_prof');
PB.labels = ["X Prior" "X Retrieved"];
PB.xlims = [2e-6,9e-6];
PB.logy = true;
PB.ydec = true;
PB.ylims = [1,1000];
%PB.ylims = [20,80];
%plot wv ppmv
plot_basic([Xb,X_ret],P, PB);

%plot T vs P
PB = PB_atm;
PB.logy = true;
PB.ydec = true;
PB.ylims = [1,1000];
PB.savestring = append(savefile,'temperature');
PB.xlab = "Temperature (K)";
PB.xlims = [50,300];
PB.labels = ["temperature"];
PB.xlims = [50,300];
plot_basic([T],P, PB);

%plot z vs P
PB.xlab = "Height (m asl)";
PB.savestring = append(savefile,"altitude");
PB.labels = ["altitude"];
PB.xlims = [10000,100000];
plot_basic([Z],P, PB);

%plot measurement response
PB = PB_atm;
PB.logy = true;
PB.ydec = true;
PB.ylims = [1,1000];
PB.savestring = append(savefile,'measurement_response');
PB.xlab = 'Measurement Response';
PB.xlims = [0,1.2];
PB.labels = ["Measurement Response"];
plot_basic([mr], P, PB);

%several averaging kernels
PB = PB_atm;
PB.logy = true;
PB.ydec = true;
PB.multi_y = true;
PB.ylims = [1,1000];
PB.savestring = append(savefile,'av_kernel');
PB.xlab = 'Averaging Kernel';
PB.xlims = [0, 0.1];
PB.labels = ["Averaging Kernel"];
av_k_plots = [];
for j = 1:5:200
    av_k_plots = [av_k_plots; A(j,:)];
end
plot_basic(av_k_plots', P, PB);


function plot_basic(x,y,PB)
    %expected fields
    pb_fields = ["xlab" "ylab" "savestring" "logy" "ydec" "labels" "xlims" "ylims" "multi_y"];
    
    %set expected fields to false if not found
    for i = 1:numel(pb_fields)
      if ~ isfield(PB, pb_fields{i})
        PB.(pb_fields{i}) = false;
      end
    end
    %initialise fig
    fig = figure('visible','off');
    hold on
    grid on

    %iterate through lines to plot
    for i = 1:length(PB.labels)
       if PB.multi_y == 1
           plot(x,y(:,i),'LineWidth',3);
       else
           plot(x(:,i),y,'LineWidth',3);
       end
    end

    %set axes variables
    if PB.logy == 1; set(gca, 'YScale', 'log'); end
    if PB.ydec == 1; set( gca, 'YDir', 'reverse' ); end
    if ~all(PB.xlims == 0); xlim(PB.xlims); end
    if ~all(PB.ylims == 0); ylim(PB.ylims); end
    %set label details
    xlabel(PB.xlab);
    ylabel(PB.ylab);
    legend(PB.labels);
    set(gca,'FontSize',16)
    set(gca,'LineWidth',.5)

    %want to add meta-data to each fig, not sure about a good way to do
    %this
    saveas(fig,append(PB.savestring,'.png'), 'png' )
    close all
 end
 
 end
