function plot_RF_struct_rast(WN_datafile,cellIDs,structu)


datarun = load_data(WN_datafile)
datarun = load_params(datarun)
datarun = load_sta(datarun);

figure;
for icell=1:length(cellIDs)
plot_rf_fit(datarun,cellIDs(icell),'fill_color',[0,0,0],'alpha',(structu(icell)-min(structu))/(max(structu)-min(structu)),'fill',true,'edge',true,'labels',true,'edge_color',[1,0,0]);
hold on;
axis image;
set(gca,'xTick',[]);set(gca,'yTick',[]);
end

end