function MyHeatMap(data,HFig,Hsbpl,ratios)
% plots HeatMap, so far gray scale; -1 black, 0 middle gray, 1 white
% data: matrix of data to plot. Should be between -1 and 1.
% HFig: figure handle
% Hsbpl: subplot handle
% ratios: how to display checkers. 1: squared checkers, 2: equal axes,
% other values: leave as it is


data=(data+1)/2;

figure(HFig)
subplot(Hsbpl)
for i=1:size(data,1)
    for j=1:size(data,2)        
        rectangle('Position',[j-1,size(data,1)-i,1,1],'FaceColor',[1 1 1]*data(i,j),'EdgeColor','none')
    end    
end
set(gca,'Visible','off')
if ratios==1
    set(gca,'DataAspectRatio',[1,1,1])
elseif ratios==2
    if size(data,1)<size(data,2)
        set(gca,'DataAspectRatio',[1,size(data,1)/size(data,2),1])
    else
        set(gca,'DataAspectRatio',[size(data,2)/size(data,1),1,1])
    end
end