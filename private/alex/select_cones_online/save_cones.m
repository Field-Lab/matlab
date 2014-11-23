function save_cones(path2save, mapName)

global datarun stim cellType myCells mean_nnd cones

height = datarun.stimulus.stixel_height * datarun.stimulus.field_height;
width = datarun.stimulus.stixel_width * datarun.stimulus.field_width;
stixel = datarun.stimulus.stixel_width;


coord=stim.coord*stixel;

R=floor(mean_nnd*stixel/2-1);% radius

tmp=gcf;


hnew=figure;
set(hnew,'Name','Select Radius','Position',[ -888   489   700   800],'ToolBar','Figure')
hplot=subplot('position',[0.1 0.1 0.7 0.7]);

hText=uicontrol('style','text','Units','normalized','position',[0.1 0.9 0.2 0.03],...
    'string','enter radius','fontsize',24);

hText=uicontrol('style','text','Units','normalized','position',[0.5 0.9 0.2 0.03],...
    'string','to finish enter 0','fontsize',16);

hEditRadius=uicontrol('style','edit','Units','normalized','position',[0.2 0.82 0.1 0.05],...
    'string',int2str(R),'fontsize',24);

flag=1;
while ~strcmp(get(hEditRadius, 'String'), '0')
    R=str2num(get(hEditRadius,'String'));
    myMap=zeros(height,width);
    for kone=1:length(coord)
        
        C = coord(kone,:);
        
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


if ~exist(path2save,'dir')
    mkdir(path2save)
end

clear info
info.name = datarun.names.rrs_prefix;
info.cellType = cellType;
info.cells = myCells;

dlmwrite([path2save mapName, '.txt'], myMap, 'delimiter', '\t', 'newline', 'pc');
save([path2save mapName '_info'],'stim', 'info', 'cones')

close(hnew)

figure(tmp);

hSavedInfo=uicontrol('style','text', 'Units', 'Normalized','position',[0.32 0.08 0.32 0.1],...
    'string',{'saved to', path2save, ['file: ' mapName '.txt'], [mapName '_info.mat']},'fontsize',16, 'fontweight', 'bold','foregroundcolor','b');







