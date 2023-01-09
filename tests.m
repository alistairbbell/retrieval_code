function res=test

a = get_apriori_h20_from_ecmwf_zimmerwald(datetime(2015,1,1,0,32,25), '/home/alistair/MIAWARA_ret/a_priori/ecmwf_2010_2015_3_9_15_21h.nc');

%disp(a(1,1,1:100))
plot_data = reshape(a,[],1);
disp('size(plot_data)')

disp(size(plot_data))
xdata = (1:1500);
plot(a(:,1), xdata )
end



