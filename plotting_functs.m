retrievalPlotting(L2, settings_metadata, savefile)
%%%%% functions to plot and save the retrieval data %%%%

%define the spectra variables to plot
freq = L2.Y.F; %frequency vector
y_sim = L2.yf; %retrieval with baseline fit
y_obs = L2.Y.Y; %original observation
bl = L2.bl; %baseline fit of spectra
res = y_obs - y_sim; %residuals

%define retrieval quantities
P = L2.species1_p;
z = L2.species1_z;
Xb = L2.species1_xa; % wv background 
X_ret = L2.species1_x;% wv retrieval
mr = L2.species1_mr; %measurement response

%create date subfolder (might change later to plot type)

plot_basic(freq/1e9, [y_sim, y_obs], 'Frequency (GHz)', 'Brightness Temperature (K)',...
    append(savefile,'Tb_y_yret'), false, false, ['Observation' 'Retrieval'])



function plot_basic(x,y,xlab, ylab, savestring, logy , ydec, labels)
    hold on
    grid on
    for i = 1:length(y(1,:))
      plot(x,y(:,i));
    end
    if logy == 1
        set(gca, 'YScale', 'log');
    end
    if ydec == 1
        set( hAxes, 'YDir', 'reverse' );
    end
    xlabel(xlab);
    ylabel(ylab);
    legend(labels);
    fn = fieldnames(settings_metadata);
    for k=1:numel(fn)
        set(gcf, fn{k} , settings_metadata.(fn{k}));
    end
    print(gcf, '-dpdf', '-append', append(savestring,'.pdf'));
 end

