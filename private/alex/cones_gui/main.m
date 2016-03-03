
clear
close all
clear
global myCells datarun cones cellType cell_indices coord_tform ctr rad fit_angle
global stim STAplotPosition RFfitsPlotPosition deleted_cells_list
cones={}; myCells=[]; deleted_cells_list=[];


run = 'data001';
path2load = ['/Volumes/Analysis/2016-01-05-1/d00-06-norefit/data001/data001'];
path2save=fullfile(server_path(), date, '/stimuli/maps_test/');
frame=4;


%% main body

mapName = ['map_', run];
cellType=4;
nnd_scale=3;
field_size = [400 400];

datarun = load_data(path2load);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all','keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

ncells = length(datarun.cell_ids);

%% preprocess RF fits
coord_tform = coordinate_transform(datarun,'sta');
ctr=nan(ncells,2);
rad=ctr;
fit_angle=nan(ncells,1);
for cell_index = 1:ncells
    the_fit = datarun.stas.fits{cell_index};
    if ~isempty(the_fit)
        ctr(cell_index,:) = the_fit.mean;
        rad(cell_index,:) = the_fit.sd;
        fit_angle(cell_index)=the_fit.angle;
    end
end


%% preprocess STA

all_sta = zeros([field_size, ncells]);
summed_rfs = 0;

brightest_star = zeros(11);
for rgc=1:ncells
    sta = datarun.stas.stas{rgc}(:,:,1,frame);
    
    if abs(min(sta(:)))>max(sta(:)) % OFF cell
        sta = -sta;
    end
    all_sta(:,:,rgc) = sta;
    t = robust_std(sta(:))*5;
    sigstix = sum(sta(:)>t);
    if sigstix>0 && sigstix<100 % possibly cell with cones
        %             sta = imresize(sta,2);
        tmp = sta;
        tmp = tmp/max(tmp(:));
        t = robust_std(tmp(:))*4;
        tmp(tmp<t) = 0;
        [r,c] = find(tmp);
        tmp1 = squareform(pdist([r c]));
        tmp1(tmp1==0)=760;
        p = find(min(tmp1)>25);
        for i=1:length(p)
            tmp(r(p(i)), c(p(i))) = 0;
        end        
        summed_rfs = summed_rfs+tmp;
        [r,c] = find(sta==max(sta(:)),1);
        if r>10 && r<field_size(1)-10 && c>10 && c<field_size(2)-10
            brightest_star = brightest_star+sta(r-5:r+5,c-5:c+5);
        end
    end
end

figure
imagesc(brightest_star)


%% local max cone finding

cone_peaks = find_local_maxima(summed_rfs, 'radius', 1, 'thresh', 0.7, 'return', 'indices');
figure
colormap gray
imagesc(summed_rfs)
hold on
plot(cone_peaks(:,2), cone_peaks(:,1), 'xg')
set(gca, 'dataaspectratio',[1 1 1])

%% GUI

hmain=figure;
set(hmain,'Name','Main Window','Position',[1601 177 1280 929],'ToolBar','Figure')

STAplotPosition=[0.5 0.46 0.45 0.45];
sta_plot=subplot('position', STAplotPosition);

current_cell = 1;

hNextCell=uicontrol('style','pushbutton','Units','normalized','position',[0.1 0.95 0.12 0.03],...
    'string','next cell','fontsize',16,'callback','show_next_sta(all_sta(current_cell+1), current_cell+1, cone_list, sta_plot)');

hAdd = uicontrol('style','pushbutton','Units','normalized','position',[0.45 0.95 0.12 0.03],...
    'string','Add','fontsize',16,'callback','add_cone(cone_list, contour_list, sta, sta_plot)');

hReject = uicontrol('style','pushbutton','Units','normalized','position',[0.6 0.95 0.12 0.03],...
    'string','Reject','fontsize',16,'callback','reject_cone()');

hAdjust

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

