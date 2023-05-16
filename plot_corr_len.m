function plot_corr_len(corr_len_2d, savefile)


%%%%% functions to plot and save the retrieval data %%%%
try
    mkdir(savefile)
catch
    'Directory already exists'
end

P = corr_len_2d(:,1);
len_scale = corr_len_2d(:,2);

%%%   define struct for atmosphere plots %%%
PB_atm.ylab =  'Pressure (Pa)';
PB_atm.logy = true;
PB_atm.ydec = true;
PB_atm.xlims = false;

%plot measurement response
PB = PB_atm;
PB.ylims = [1,1000];
PB.savestring = append(savefile,'correlation_length');
PB.xlab = 'Correlation Length scale';
PB.labels = ["Correlation Length"];
plot_basic([len_scale], P, PB);

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