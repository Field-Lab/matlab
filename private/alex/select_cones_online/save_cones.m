function save_cones(path2save, mapName)

global datarun stim cellType myCells mean_nnd cones myMap conesSavePosition
global myCluster
myCluster=[];

height = datarun.stimulus.stixel_height * datarun.stimulus.field_height;
width = datarun.stimulus.stixel_width * datarun.stimulus.field_width;
stixel = datarun.stimulus.stixel_width;


coord=stim.coord*stixel;

R=max(1, floor(mean_nnd*stixel/2-1));% radius

tmp=gcf;

conesSavePosition=[0.1 0.1 0.7 0.7];
hnew=figure;
set(hnew,'Name','Select Radius','Position',[ -888   489   700   800],'ToolBar','Figure')
hplot=subplot('position',[0.1 0.1 0.7 0.7]);

hText=uicontrol('style','text','Units','normalized','position',[0.1 0.9 0.25 0.07],...
    'string',{'enter radius', 'to exit enter 0'},'fontsize',24);

% hText=uicontrol('style','text','Units','normalized','position',[0.1 0.9 0.2 0.03],...
%     'string','to finish enter 0','fontsize',16);

hEditRadius=uicontrol('style','edit','Units','normalized','position',[0.2 0.82 0.1 0.05],...
    'string',int2str(R),'fontsize',24);

% select radius of stim flash
flag=1;
while ~strcmp(get(hEditRadius, 'String'), '0')
    R=str2num(get(hEditRadius,'String'));
    myMap=zeros(height,width);
    for kone=1:length(coord)
        
        C = coord(kone,:);
        if R==1 % if single pixel stimulus desired
            myMap(round(C(1)),round(C(2)))=kone;            
        else            
            t = linspace(0, 2*pi, 100);
            x = round(R*cos(t) + C(1));
            y = round(R*sin(t) + C(2));
            %     figure
            %     plot(x,y)
            cc=zeros((max(x)-min(x)+1)*(max(y)-min(y)+1),2);
            cnt=1;
            for i=min(x):max(x)
                for j=min(y):max(y)
                    cc(cnt,:)=[i,j];
                    cnt=cnt+1;
                end
            end
            IN = inpolygon(cc(:,1),cc(:,2),x,y);
            myMask=cc(IN,:);
            myMask(myMask<1)=1;
            myMask(myMask(:,1)>size(myMap,1),1)=size(myMap,1);
            myMask(myMask(:,2)>size(myMap,2),2)=size(myMap,2);
            
            for i=1:length(myMask)
                myMap(myMask(i,1),myMask(i,2))=kone;
            end
        end

        
    end
    myMap=myMap';
    
    xx=get(gca, 'XLim');
    yy=get(gca, 'YLim');
    imagesc(myMap)
    set(gca, 'DataAspectRatio',[1 1 1]);
    if flag
        axis tight
        flag=0;
    else
        axis([xx yy]);
    end
    
    waitfor(hEditRadius,'string')
end

stim.radius=R;
stim.stixel=stixel;


% select mode of stimulation
hTemplate = uicontrol('Style','pushbutton','Units','normalized','Position',[0.6 0.92 0.15 0.05],...
    'String','Do clusters','fontsize',16, 'callback','select_clusters(cones)');
hTemplate1 = uicontrol('Style','checkbox','Units','normalized','Position',[0.6 0.85 0.15 0.05],...
    'String','Independent','visible','on');

hSave=uicontrol('style','togglebutton','Units','normalized','position',[0.8 0.85 0.15 0.05],...
    'string','save','fontsize',16);

waitfor(hSave,'Value',1)

if ~isempty(myCluster)
    i=1;
    while ~isempty(myCluster{i})
        i=i+1;
    end
    if i>1
        myCluster=myCluster(1:i-1);
        
        if get(hTemplate1,'Value') % independent stimulation
            
            for i=1:length(myCluster)
                for j=1:length(myCluster{i})
                    myMap(myMap==myCluster{i}(j))=-j;
                end
            end
            
        else % all cones in cluster flash simultaneously
            for i=1:length(myCluster)
                for j=1:length(myCluster{i})
                    myMap(myMap==myCluster{i}(j))=-i;
                end
            end
        end
        
        % delete unclustered cones
        myMap(myMap>0)=0;
        
        % make clustered cones positive
        myMap=-myMap;
        myMap(abs(myMap)<0.1)=0;
    end
end

if ~exist(path2save,'dir')
    mkdir(path2save)
end

info.name = datarun.names.rrs_prefix;
info.cellType = cellType;
info.cells = myCells;

file_path=[path2save mapName '.txt'];

dlmwrite(file_path, myMap, 'delimiter', '\t', 'newline', 'pc');
save([path2save mapName '_info'],'stim', 'info', 'cones')

close(hnew)

figure(tmp);

% hSavedInfo=uicontrol('style','text', 'Units', 'Normalized','position',[0.32 0.08 0.32 0.1],...
%     'string',{'saved to', path2save, ['file: ' mapName '.txt'], [mapName '_info.mat']},'fontsize',16, 'fontweight', 'bold','foregroundcolor','b');


savedMap=dlmread(file_path);
figure
imagesc(savedMap)
title(file_path,'Interpreter','None','fontsize',16, 'fontweight', 'bold')
disp(file_path)











