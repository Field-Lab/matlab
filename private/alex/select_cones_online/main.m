% load data
clear
close all
clear
global myCells datarun cones cellType cell_indices coord_tform ctr rad fit_angle
global stim STAplotPosition RFfitsPlotPosition
cones={}; myCells=[];

date='2013-08-19-2';
path2load = fullfile(server_path(), date, '/streamed/data001/data001');
path2save=fullfile(server_path(), date, '/stimuli/maps/');

date='2011-10-25-5';
run='data001-0';
path2load = fullfile(server_path(), [date, '/streamed/',run,'/',run]);
path2save=fullfile(server_path(), date, '/stimuli/maps/');

date='2011-12-13-2';
run='data008-0';
path2load = fullfile(server_path(), [date, '/streamed/',run,'/',run]);
path2save=fullfile(server_path(), date, '/stimuli/maps/');


date = '2012-09-13-2';
run = 'data009';
path2load = fullfile(server_path(), [date, '/',run,'/',run]);
path2save=fullfile(server_path(), date, '/stimuli/maps/');
frame=5;

date = '2013-10-10-5';
run = 'data002';
path2load = fullfile(server_path(), [date, '/streamed/',run,'/',run]);
path2save=fullfile(server_path(), date, '/stimuli/maps/');
frame=6;

date = '2010-03-05-2';
run = 'data001';
path2load = fullfile(server_path(), [date, '/',run,'/',run]);
path2save=fullfile(server_path(), date, '/stimuli/maps/');
frame=5;

date = '2013-08-19-4';
run = 'data001';
path2load = fullfile(server_path(), [date, '/streamed/',run,'/',run]);
path2save=fullfile(server_path(), date, '/stimuli/maps/');
frame=5;



date = '2013-08-19-5';
run = 'data001';
path2load = fullfile(server_path(), [date, '/streamed/',run,'/',run]);
path2save=fullfile(server_path(), date, '/stimuli/maps/');
frame=5;


mapName = ['map_', run];
cellType=5;
nnd_scale=3;



datarun = load_data(path2load);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);


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

RFfitsPlotPosition=[0.01 0.46 0.45 0.45];

hplot=subplot('position',RFfitsPlotPosition);
plot_rf_summaries(datarun, {cellType}, 'clear', false, 'label', true, 'label_color', 'k', 'plot_fits', true, 'fit_color', 'k')

STAplotPosition=[0.5 0.46 0.45 0.45];
% buttons

hSelectCells=uicontrol('style','pushbutton','Units','normalized','position',[0.1 0.95 0.12 0.03],...
    'string','select cells','fontsize',16,'callback','select_cells(hplot, nnd_scale, frame)');


hZoomOut=uicontrol('style','pushbutton','Units','normalized','position',[0.30 0.95 0.08 0.03],...
    'string','zoom out','fontsize',16,'callback','show_sta(0)');

hAddManual=uicontrol('style','pushbutton','Units','normalized','position',[0.45 0.95 0.12 0.03],...
    'string','add manually','fontsize',16,'callback','handle_cones(1,length(myCells))');


hDeleteCone=uicontrol('style','pushbutton','Units','normalized','position',[0.6 0.95 0.12 0.03],...
    'string','delete cone','fontsize',16,'callback','handle_cones(2,length(myCells))');

hDeleteCell=uicontrol('style','pushbutton','Units','normalized','position',[0.75 0.95 0.12 0.03],...
    'string','delete cell','fontsize',16,'callback','delete_cell(1)');

hDeleteAllCones=uicontrol('style','pushbutton','Units','normalized','position',[0.89 0.95 0.1 0.03],...
    'string','delete all cones','fontsize',16,'callback','handle_cones(4,length(myCells))');




hAdjust=uicontrol('style','pushbutton','Units','normalized','position',[0.1 0.4 0.12 0.03],...
    'string','adjust','fontsize',16,'callback','adjust_cones(frame, hTemplate)');

hSubunits=uicontrol('style','pushbutton','Units','normalized','position',[0.40 0.38 0.09 0.03],...
    'string','subunits','fontsize',16,'callback','subunits(datarun, path2load, date)');

hCheckOtherCells=uicontrol('style','pushbutton','Units','normalized','position',[0.4 0.34 0.09 0.03],...
    'string','check others','fontsize',16,'callback','check_other_cells(frame)');

hSave=uicontrol('style','pushbutton','Units','normalized','position',[0.40 0.3 0.09 0.03],...
    'string','save','fontsize',16,'callback','save_cones(path2save, mapName)');


hTemplate(1) = uibuttongroup('visible','on','Position',[0.32 0.3 0.07 0.1]);
hTemplate(2) = uicontrol('Style','radiobutton','String','Gauss',...
    'pos',[10 10 70 30],'parent',hTemplate(1),'HandleVisibility','on');
hTemplate(3) = uicontrol('Style','radiobutton','String','Template',...
    'pos',[10 50 70 30],'parent',hTemplate(1),'HandleVisibility','on');
set(hTemplate(1),'SelectedObject',hTemplate(2)); 
set(hTemplate(1),'Visible','on');







