% load data
clear
close all
clear
global myCells datarun cones cellType cell_indices coord_tform ctr rad fit_angle

cones={}; myCells=[];

datarun = load_data(fullfile(server_path(), '2013-08-19-2/streamed/data001/data001'));

cellType=4;

datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);

nnd_scale=3;
frame=5;


cell_indices = get_cell_indices(datarun, {cellType});
coord_tform = coordinate_transform(datarun,'sta');

ctr=nan(length(datarun.cell_ids),2);
rad=ctr;
fit_angle=nan(length(datarun.cell_ids),1);
for cell_index = 1:length(datarun.cell_ids)
    % get the fit
    the_fit = datarun.stas.fits{cell_index};
    % skip if doesn't exist
    if isempty(the_fit);continue;end
    % get center
    ctr(cell_index,:) = the_fit.mean;
    % get center radius
    rad(cell_index,:) = the_fit.sd;
    fit_angle(cell_index)=the_fit.angle;
end


hmain=figure;
set(hmain,'Name','Main Window','Position',[1601 177 1280 929],'ToolBar','Figure')

hplot=subplot('position',[0.01 0.62 0.3 0.3]);
plot_rf_summaries(datarun, {cellType}, 'clear', false, 'label', true, 'label_color', 'k', 'plot_fits', true, 'fit_color', 'k')


% buttons

hCells=uicontrol('style','pushbutton','Units','normalized','position',[0.1 0.945 0.12 0.03],...
    'string','select cells','fontsize',16,'callback','select_cells(hplot, nnd_scale, frame)');


hZoomOut=uicontrol('style','pushbutton','Units','normalized','position',[0.30 0.93 0.08 0.03],...
    'string','zoom out','fontsize',16,'callback','show_sta(0)');

hAddManual=uicontrol('style','pushbutton','Units','normalized','position',[0.45 0.945 0.12 0.03],...
    'string','add manually','fontsize',16,'callback','handle_cones(1,length(myCells))');


hDeleteCone=uicontrol('style','pushbutton','Units','normalized','position',[0.6 0.945 0.12 0.03],...
    'string','delete cone','fontsize',16,'callback','handle_cones(2,length(myCells))');


hDeleteCell=uicontrol('style','pushbutton','Units','normalized','position',[0.4 0.5 0.12 0.03],...
    'string','delete cell','fontsize',16,'callback','delete_cell');


hAdjust=uicontrol('style','pushbutton','Units','normalized','position',[0.1 0.52 0.12 0.03],...
    'string','adjust','fontsize',16,'callback','adjust_cones(frame)');


%% still do!

hCheckOtherCells=uicontrol('style','pushbutton','Units','normalized','position',[0.67 0.52 0.12 0.03],...
    'string','check others','fontsize',16,'callback','check_other_cells(frame)');

hSave=uicontrol('style','pushbutton','Units','normalized','position',[0.40 0.3 0.12 0.03],...
    'string','save','fontsize',16,'callback','save_stim');





