
for letterCode='ABCDEFGH'
    letterCode
    load(['/Users/alexth/Desktop/dataruns/datarun',letterCode])
end


global myCell selected_cones selected_params

datarun=datarunF;

datarun=load_sta(datarun,'load_sta','all');

masks = datarun.cones.mosaic.voronoi_masks;
excludes = [];
indexes = 1:length(masks);
indexes = setdiff(indexes, excludes);
min_neighbor_dist = 2;
max_self_dist     = 3.5;
spaced = space_voronoi_masks(datarun, min_neighbor_dist, max_self_dist);
cone_map = index_masks(spaced, num2cell(indexes));

[v,c] = voronoin(datarun.cones.centers);


%%   
% plot ranges for sigma

titles={'ON parasol','OFF parasol','ON midget','OFF midget'};

figure
for i=1:4
    lowS=zeros(size(datarun.cell_types{i}.cell_ids));
    highS=zeros(size(datarun.cell_types{i}.cell_ids));
    cnt=1;
    for currentCell=datarun.cell_types{i}.cell_ids    
        mosaic_weights=datarun.cones.weights(:,datarun.cell_ids==currentCell);
        a=optThreshold(datarun.cell_ids==currentCell);
        if a>0
            tmp=sort(mosaic_weights,'descend');
            a=tmp(a-1);
            mosaic_weights(mosaic_weights<a)=0;
            a=a/max(mosaic_weights);
            mosaic_weights=mosaic_weights/max(mosaic_weights);
            m=find(mosaic_weights);
            
            k=[];
            sigmaWeight=1.1:0.1:10;
            for threshold=sigmaWeight
                
                [mosaic_weights, selection, extras] = select_cone_weights(datarun, currentCell,...
                    'thresh', threshold, 'radius', [0 inf], 'polarity', 1,'contiguity', false);
                m1=find(selection);
                k=[k length(setxor(m,m1))];
            end
            kk=min(k);
            k=find(k==kk);
            
            lowS(cnt)=sigmaWeight(k(1));
            highS(cnt)=sigmaWeight(k(end));
        end
        cnt=cnt+1;
    end
    subplot(2,2,i)
    plot(highS,'*-')
    hold on 
    plot(lowS,'-*')
    title(titles{i})
end

 
%plot ranges for weight

titles={'ON parasol','OFF parasol','ON midget','OFF midget'};

figure
for i=1:4
    incl=zeros(size(datarun.cell_types{i}.cell_ids));
    excl=zeros(size(datarun.cell_types{i}.cell_ids));
    cnt=1;
    for currentCell=datarun.cell_types{i}.cell_ids    
        mosaic_weights=datarun.cones.weights(:,datarun.cell_ids==currentCell);
        a=optThreshold(datarun.cell_ids==currentCell);
        if a>0
            tmp=sort(mosaic_weights,'descend');
            tmp=tmp/max(mosaic_weights);
            incl(cnt)=tmp(a-1);
            excl(cnt)=tmp(a);
        end
        cnt=cnt+1;
    end
    subplot(2,2,i)
    plot(incl,'*-')
    hold on 
    plot(excl,'-*')
    title(titles{i})
end


 
%plot ranges for weight - choose datarun

datarun=datarunC;

titles={'ON parasol','OFF parasol','ON midget','OFF midget'};

figure
set(gcf,'Name',datarun.names.short_name)
for i=1:4

    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
    subplot(2,2,i)
    plot(datarun.LowestInclWeight(selectRGCind),'*-')
%     plot(datarun.LowestInclWeight(selectRGCind)./datarun.RobustSTD(selectRGCind),'*-')
    hold on 
    plot(datarun.HighestExclWeight(selectRGCind),'-*')
%     plot(datarun.HighestExclWeight(selectRGCind)./datarun.RobustSTD(selectRGCind),'-*')
    title(titles{i})
end


datarun=datarunD;

titles={'ON parasol','OFF parasol','ON midget','OFF midget'};

figure
set(gcf,'Name',datarun.names.short_name)
for i=1:4

    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
    subplot(2,2,i)
    plot(datarun.RobustSTD(selectRGCind),'*-')
%     plot(datarun.LowestInclWeight(selectRGCind)./datarun.RobustSTD(selectRGCind),'*-')
    hold on 
    plot(datarun.STD(selectRGCind),'-*r')
%     plot(datarun.HighestExclWeight(selectRGCind)./datarun.RobustSTD(selectRGCind),'-*')
    title(titles{i})
end

datarun=datarunC;

titles={'ON parasol','OFF parasol','ON midget','OFF midget'};

figure
set(gcf,'Name',datarun.names.short_name)
for i=1:4

    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
    subplot(2,2,i)
    plot(datarun.LowestInclWeight(selectRGCind)./datarun.STD(selectRGCind),'*-')
%     plot(datarun.LowestInclWeight(selectRGCind)./datarun.RobustSTD(selectRGCind),'*-')
    hold on 
%     plot(datarun.STD(selectRGCind),'-*r')
%     plot(datarun.HighestExclWeight(selectRGCind)./datarun.RobustSTD(selectRGCind),'-*')
    title(titles{i})
end



%%

myFig=figure;
set(myFig,'Name',datarun.names.short_name,'Units','normalized','Position',[1.0000 0 0.8000 0.7925])

handles.cellType=uibuttongroup('visible','on','Units','normalized','Position',[0.01 0.7 .1 0.2]);
handles.onP=uicontrol('style','Radio','String','ON parasol','FontWeight','bold','Units','normalized','position',[0.02 0.7 0.9 0.2],...
    'parent',handles.cellType,'HandleVisibility','on','Value',1);
handles.offP=uicontrol('style','Radio','String','OFF parasol','FontWeight','bold','Units','normalized','position',[0.02 0.5 0.9 0.2],...
    'parent',handles.cellType,'HandleVisibility','on','Value',0);
handles.onM=uicontrol('style','Radio','String','ON midget','FontWeight','bold','Units','normalized','position',[0.02 0.3 0.9 0.2],...
    'parent',handles.cellType,'HandleVisibility','on','Value',0);
handles.offM=uicontrol('style','Radio','String','OFF midget','FontWeight','bold','Units','normalized','position',[0.02 0.1 0.9 0.2],...
    'parent',handles.cellType,'HandleVisibility','on','Value',0);
set(handles.cellType,'SelectionChangeFcn','plot_new(datarun,datarun.cell_types{5-find(cell2mat(get(get(handles.cellType,''children''),''value'')))}.cell_ids,v,c,1,get(h8,''string''))');


h7=uicontrol('style','pushbutton','Units','normalized','position',[0.01 0.61 0.1 0.07],'string','show STA',...
    'callback','show_sta_IR(datarun,v,c,myCell)');

h4=uicontrol('style','togglebutton','Units','normalized','position',[0.01 0.5 0.1 0.07],'string','auto zoom',...
    'Value',1,'callback','auto_zoom(datarun,myCell,get(h4,''value''),get(h5,''string''))');

h5=uicontrol('style','edit','Units','normalized','position',[0.06 0.57 0.04 0.03],'string','3','fontsize',16);
h6=uicontrol('style','text','Units','normalized','position',[0.01 0.57 0.04 0.03],'string','zoom, sigmas','fontsize',12);

h3=uicontrol('style','togglebutton','Units','normalized','position',[0.01 0.4 0.1 0.07],'string','contig True',...
    'Value',1,'callback','recalc_cone_selection(datarun,myCell,get(h,''Value''),v,c,get(h3,''Value''))');

h2=uicontrol('style','togglebutton','Units','normalized','position',[0.01 0.3 0.1 0.07],'string','voronoi ON',...
    'Value',1,'callback','switch_voronoi(v,c,get(h2,''Value''))');

h1=uicontrol('style','pushbutton','Units','normalized','position',[0.01 0.2 0.1 0.07],'string','next cell',...
    'callback','plot_new(datarun,datarun.cell_types{5-find(cell2mat(get(get(handles.cellType,''children''),''value'')))}.cell_ids,v,c,0,get(h8,''string''))');

h=uicontrol('style','slider','Units','normalized','position',[0.1 -0.15 0.7 0.2],'min',0.05,'max',15.05,'value',10,...
    'SliderStep',[1/300,1/15],'callback','recalc_cone_selection(datarun,myCell,get(h,''Value''),v,c,get(h3,''Value''))');
 
h8=uicontrol('style','edit','Units','normalized','position',[0.06 0.15 0.04 0.03],'string','5','fontsize',16,...
'callback','plot_new(datarun,datarun.cell_types{5-find(cell2mat(get(get(handles.cellType,''children''),''value'')))}.cell_ids,v,c,-1,get(h8,''string''))');
h9=uicontrol('style','text','Units','normalized','position',[0.01 0.15 0.04 0.03],'string','frame','fontsize',12);

h10=uicontrol('style','pushbutton','Units','normalized','position',[0.85 0.02 0.07 0.07],'string','clear marks',...
    'callback','recalc_cone_selection(datarun,myCell,get(h,''Value''),v,c,get(h3,''Value''),true)');




s=uicontrol('style','pushbutton','Units','normalized','position',[0.2 0.95 0.07 0.04],'string','save opt',...
    'fontweight','bold','fontsize',12,'callback','save_results(datarun,myCell,1,get(h,''Value''),get(h3,''Value''))');

s2=uicontrol('style','pushbutton','Units','normalized','position',[0.35 0.95 0.07 0.04],'string','save min',...
    'fontweight','bold','fontsize',12,'callback','save_results(datarun,myCell,2,get(h,''Value''),get(h3,''Value''))');

s3=uicontrol('style','pushbutton','Units','normalized','position',[0.5 0.95 0.07 0.04],'string','save max',...
    'fontweight','bold','fontsize',12,'callback','save_results(datarun,myCell,3,get(h,''Value''),get(h3,''Value''))');
