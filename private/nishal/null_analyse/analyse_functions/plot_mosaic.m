
function plot_mosaic(datarun,InterestingCell_vis_id,ref_cell_number)
h=figure('Color','w');


cellID_use=InterestingCell_vis_id(ref_cell_number);
plot_rf_fit(datarun, 'On Parasol','fill',false,'edge',true,'labels',false,'edge_color',[0,0,1]);
hold on
plot_rf_fit(datarun, 'Off Parasol','fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
hold on
%plot_rf_fit(datarun, cellID_use,'fill_color',[0,1,0],'alpha',0.3,'fill',true,'edge',false,'labels',true,'edge_color',[1,0,0]);
plot_rf_fit(datarun,cellID_use,'fill_color',[0,1,0],'alpha',0.3,'fill',true,'edge',false,'labels',false,'edge_color',[1,0,0]);

%xlim([0,Filtdim1]);
%ylim([0,Filtdim2]);

title('Moasic (Blue: On Parasol, Red: Off Parasol)')
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'YColor','w')
set(gca,'XColor','w')
end