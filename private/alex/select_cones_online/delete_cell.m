function delete_cell

global myCells cones myLimits ctr rad fit_angle coord_tform datarun

% delete all cones form the STA plot
hPlot=findobj('position',[0.4 0.62 0.3 0.3]);
subplot(hPlot)
for i=1:length(cones{end})
    children = get(gca, 'children');
    myPoint=findobj('XData',cones{end}(i,1),'YData',cones{end}(i,2), 'color','r');
    delete(children(children==myPoint));
end

% set new axis
axis(myLimits)

% recolor rf fit
hPlot=findobj('position',[0.01 0.62 0.3 0.3]);
subplot(hPlot)
datInd=find(datarun.cell_ids==myCells(end));
[X, Y] = drawEllipse([ctr(datInd,:) rad(datInd,:) fit_angle(datInd)]);
[X, Y] = tformfwd(coord_tform, X, Y);
hold on
plot(X,Y,'k')

% update this cell info
uicontrol('style','text', 'Units', 'Normalized','position',[0.75 0.8 0.2 0.15],...
    'string',{'','',['Cell ',int2str(myCells(end)), ' was deleted!']},...
    'BackgroundColor',[.7 .9 .7],'ForegroundColor','r','fontsize',24, ...
    'fontweight', 'bold');

% delete from arrays
myCells(end)=[];
cones(end)=[];

% update all cones plot
show_all_cones;

% update cone info
update_cone_info;


