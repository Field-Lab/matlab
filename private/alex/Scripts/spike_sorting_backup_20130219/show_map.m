function show_map
global rgc uicontr_units spikes pc1 pc2


%check if at least 1 unit is selected
% check if rgc is not empty
col=[1,0,0;0,1,0;0.5,0.5,0;0.25,0.75,0.75;0.75,0.25,0.75;0.05,0.58,0.05;0.9,0.45,0.1;0.4,0.25,0.75;0.86,0.86,0];

selected_units=[];
for j=1:size(rgc,2)
    if get(uicontr_units(j),'Value')==1
        selected_units=[selected_units j];
    end
end

figure;
set(gcf,'Name',['Feature map of units ',int2str(j)],'Units','normalized','Position',[0.04 0.03 0.81 0.85])
axes('DrawMode','fast');

for i=1:15
    subplot(3,5,i)
    if i==1
        xdata=pc1;
        ydata=1:length(pc1);
        xdataName='pc1';
        ydataName='time';
    elseif i==2
        xdata=pc2;
        ydata=1:length(pc1);        
        xdataName='pc2';
        ydataName='time';
    elseif i==3
        xdata=spikes(13,:);
        ydata=1:length(pc1);        
        xdataName='13';
        ydataName='time';
    elseif i==4
        xdata=pc1;
        ydata=pc2;        
        xdataName='pc1';
        ydataName='pc2';
    elseif i==5
        xdata=min(spikes);
        ydata=max(spikes);        
        xdataName='min';
        ydataName='max';
    elseif i==6
        xdata=pc1;
        ydata=min(spikes);          
        xdataName='pc1';
        ydataName='min';
    elseif i==7
        xdata=pc2;
        ydata=min(spikes);
        xdataName='pc2';
        ydataName='min';
    elseif i==8
        xdata=pc1;
        ydata=max(spikes);
        xdataName='pc1';
        ydataName='max';
    elseif i==9
        xdata=pc2;
        ydata=max(spikes);
        xdataName='pc2';
        ydataName='max';
    elseif i==10
        xdata=spikes(6,:);
        ydata=spikes(13,:);
        xdataName='6';
        ydataName='13';
    elseif i==11
        xdata=spikes(8,:);
        ydata=spikes(13,:);
        xdataName='8';
        ydataName='13';
    elseif i==12
        xdata=spikes(17,:);
        ydata=spikes(13,:);
        xdataName='17';
        ydataName='13';
    elseif i==13
%         xdata=spikes(19,:);
%         ydata=spikes(13,:);
%         xdataName='19';
%         ydataName='13';
        xdata=max(spikes)./min(spikes);
        ydata=1:length(pc1);
        xdataName='max/min';
        ydataName='time';
    elseif i==14
        xdata=spikes(21,:);
        ydata=spikes(13,:);
        xdataName='21';
        ydataName='13';
    elseif i==15
        xdata=spikes(8,:);
        ydata=spikes(21,:);
        xdataName='8';
        ydataName='21';
    end

        
    hold on

    for k=1:length(selected_units)
        plot(xdata(rgc{selected_units(k)}),ydata(rgc{selected_units(k)}),'.','MarkerSize',0.1,'color',col(selected_units(k),:))
    end
    title([xdataName,' vs ', ydataName],'Fontsize',6)
    set(gca,'XTickLabel','','YTickLabel','')
    axis tight
end



end


% default map:
% pc1 vs time
% pc2 vs time