%% load data

date = '2015-10-6-5';
run = 'data006';
path2load = '/Volumes/Acquisition/Analysis/2015-10-06-5/data006/data006';
path2save=fullfile(server_path(), date, '/stimuli/maps/');
frame=6;


%%


datarun = load_data(path2load);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);
celllist = datarun.cell_ids;

%%

global acceptedCones

hmain=figure;
set(hmain,'Name','Main Window','Position',[-1702 146 1280 929],'ToolBar','Figure')

ZoomSTAplotPosition=[0.55 0.55 0.42 0.42];
hZoomSTA=subplot('position',ZoomSTAplotPosition, 'DataAspectRatio', [1 1 1]);

STAplotPosition=[0.55 0.05 0.42 0.45];
hSTA=subplot('position',STAplotPosition, 'DataAspectRatio', [1 1 1]);

hText=uicontrol('style','text','Units','normalized','position',[0.32 0.87 0.06 0.03],...
    'string',{'radius'},'fontsize',24);

hRadius=uicontrol('style','edit','Units','normalized','position',[0.32 0.84 0.06 0.03],...
    'string','1','fontsize',24);

hCellList=uicontrol('style','listbox','Units','normalized','position',[0.4 0.1 0.08 0.85],...
    'string',{int2str(celllist')},'fontsize',16,'callback', 'newc = display_sta(get(hCellList, ''Value''), datarun, hSTA, hZoomSTA, get(hRadius, ''String''))');

hDeleteCell=uicontrol('style','pushbutton','Units','normalized','position',[0.75 0.95 0.12 0.03],...
    'string','delete cell','fontsize',16,'callback','delete_cell(1)');

