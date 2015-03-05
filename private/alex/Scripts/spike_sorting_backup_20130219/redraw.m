function redraw(xID,yID)

global spikes pc1 pc2 tdata unitPoly rgc isi_all uicontr_units src...
    getSlice1 getSlice2 visibleUnits fatPoints signUnit otherSignUI...
    getFatUI timeUnit ifUpdateISI keepAxesUI ifUndo subUnit ...
    distance2template getTemplateX getTemplateY selectValidUI
persistent sbpl_cnt 

figure(src)
switch xID
    case 1
        xdata=pc1;
        xlab='pc1';
    case 2
        xdata=pc2;
        xlab='pc2';
    case 3
        a=get(getSlice1,'String');
        xdata=spikes(str2num(a),:);
        xlab=['slice 1: point ',a];
    case 4
        a=get(getSlice2,'String');
        xdata=spikes(str2num(a),:);
        xlab=['slice 2: point ',a];
    case 5
        xdata=tdata;
        xlab='time (spike N)';
    case 6
        xdata=min(spikes);
        xlab='minima';
    case 7
        xdata=max(spikes);
        xlab='maxima';
    case 8 % template X
        a=get(getTemplateX,'String');
        xdata=distance2template(:,str2num(a));
        xlab=['template ',a];
end
switch yID
    case 1
        ydata=pc1;
        ylab='pc1';
    case 2
        ydata=pc2;
        ylab='pc2';
    case 3
        a=get(getSlice1,'String');
        ydata=spikes(str2num(a),:);
        ylab=['slice 1: point ',a];
    case 4
        a=get(getSlice2,'String');
        ydata=spikes(str2num(a),:);
        ylab=['slice 2: point ',a];        
    case 5
        ydata=tdata;
        ylab='time (spike N)';
    case 6
        ydata=min(spikes);
        ylab='minima';
    case 7
        ydata=max(spikes);
        ylab='maxima';
    case 8 % template Y
        a=get(getTemplateY,'String');
        ydata=distance2template(:,str2num(a));
        ylab=['template ',a];
end

unitsTraced=size(rgc,2)*(~isempty(rgc{1}));

col=[1,0,0;0,1,0;0.5,0.5,0;0.25,0.75,0.75;0.75,0.25,0.75;0.05,0.58,0.05;0.9,0.45,0.1;0.4,0.25,0.75;0.86,0.86,0];

axesRange=[get(gca,'XLim') get(gca,'YLim')];

if ~get(selectValidUI,'Value')
    plot(xdata,ydata,'.','MarkerSize',0.1)
    hold on
end

for k=1:unitsTraced
    while length(visibleUnits)<k
        visibleUnits(end+1)=1;
    end
    if visibleUnits(k)==1
        if signUnit(k)
            markerType=get(otherSignUI,'String');
            if markerType==''
                markerType='.';
            end
        else
            markerType='.';
        end
        if fatPoints(k)
            fatness=str2num(get(getFatUI,'String'));
        else
            fatness=0.1;
        end
        if ~get(selectValidUI,'Value') || ~isempty(intersect(subUnit,rgc{k}))
            plot(xdata(rgc{k}),ydata(rgc{k}),markerType,'MarkerSize',fatness,'color',col(mod(k-1,9)+1,:));
            hold on
        end
    end
end

% if new unit - plot
if size(unitPoly,1)>2
    IN = insidepoly(xdata,ydata,unitPoly(:,1),unitPoly(:,2));
    rgc{unitsTraced+1}=tdata(IN==1)';
    if get(selectValidUI,'Value')
        rgc{unitsTraced+1}=intersect(subUnit,rgc{unitsTraced+1});
    end
    plot(xdata(rgc{unitsTraced+1}),ydata(rgc{unitsTraced+1}),'.','MarkerSize',0.1,'color',col(mod(unitsTraced,9)+1,:))
    unitPoly=[];
    visibleUnits(unitsTraced+1)=1;
    fatPoints(unitsTraced+1)=0;
    signUnit(unitsTraced+1)=0;
    ifUpdateISI=1;    
elseif ~isempty(unitPoly)
    display('Polygon should have at least 3 points!')
end
hold off
if get(keepAxesUI,'Value')==1
    axis(axesRange)
end
xlabel(xlab)
ylabel(ylab)


if ifUpdateISI
    if ~ifUndo
        undo(0,1,1);
    else
        ifUndo=0;
    end
    
    % show isi
    
    h = findobj('Name','ISI <5ms; % - less than 2ms');
    if ~isempty(h)
        for i=1:length(sbpl_cnt)
            delete(sbpl_cnt(i))            
        end
        clear sbpl_cnt        
    elseif ~isempty(rgc{1})
        h=figure('Name','ISI <5ms; % - less than 2ms','Units','normalized','position',[0.206 0.5171 0.307 0.4075]);
    end
    
    if ~isempty(rgc{1})
        figure(h)
        if length(rgc)>9
            display('You really, really wanna more than 9 units?... code a new figure ^^')
        end
     
        n=ceil(length(rgc)/3);
        if length(rgc)==1
            l=1;
        elseif length(rgc)==2 || length(rgc)==4
            l=2;
        else
            l=3;
        end
        
        isi=cell(length(rgc),1);
        for k=1:min(length(rgc),9)
            isi{k}=diff(1000/str2num(timeUnit)*isi_all(rgc{k})); % in ms
            sbpl_cnt(k)=subplot(n,l,k);
            hist(isi{k}(isi{k}<5),10)
            h = findobj(gca,'Type','patch');
            set(h,'FaceColor',col(mod(k-1,9)+1,:),'EdgeColor','k')
            a=num2str(sum(isi{k}<=2)/numel(isi{k})*100+0.001);
            set(gca,'ytick',0,'yticklabel','')
            title([a(1:4),'% N=',int2str(numel(isi{k}))],'FontSize',10)
        end
    else
        close(h);
    end
    
    % show units mean and std
    if ~isempty(rgc{1})
        show_spikes(0)
    else
        t = findobj('Name','Units mean');
        if ~isempty(t)
            close(t);
        end
    end    
end

ifUpdateISI=0;

% manage unit selection uicontrols
figure(src)
if ~isempty(rgc{1}) && sum(size(uicontr_units,1)~=size(rgc,2)) % at least one rgc; amount of rgc and uicontrols is not equal
    for i=numel(uicontr_units):-1:1
        delete(uicontr_units(i));        
    end
    uicontr_units=zeros(length(rgc),2);
    for i=1:length(rgc)
        uicontr_units(i,1)=uicontrol('style','checkbox');        
        set(uicontr_units(i,1),'Units','normalized','position',[0.1+(i-1)*0.08 0.95 0.02 0.016],'BackgroundColor',[0.8 0.8 0.8])
        %'position',[50*i 715 15 15]);
        uicontr_units(i,2)=uicontrol('style','text','String',['Unit ',int2str(i)]);
        set(uicontr_units(i,2),'Units','normalized','position',[0.12+(i-1)*0.08 0.948 0.04 0.018],'Backgroundcolor',col(mod(i-1,9)+1,:));
    end
elseif isempty(rgc{1}) && ~isempty(uicontr_units) % all rgc have been deleted
    for i=numel(uicontr_units):-1:1
        delete(uicontr_units(i));        
    end
    uicontr_units=[];
end

set(src,'pointer','arrow')