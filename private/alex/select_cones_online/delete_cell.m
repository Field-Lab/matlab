function delete_cell(flag)

global myCells cones myLimits ctr rad fit_angle coord_tform datarun
global STAplotPosition RFfitsPlotPosition
persistent deleted_cells_list

if isempty(deleted_cells_list)
    deleted_cells_list=[];
end

if flag
    
    % delete all cones form the STA plot
    hPlot=findobj('position',STAplotPosition);
    subplot(hPlot)
    for i=1:length(cones{end})
        children = get(gca, 'children');
        myPoint=findobj('XData',cones{end}(i,1),'YData',cones{end}(i,2), 'color','r');
        delete(children(children==myPoint));
    end
    
    % set new axis
    axis(myLimits)
    
    % recolor rf fit
    hPlot=findobj('position',RFfitsPlotPosition);
    subplot(hPlot)
    datInd=find(datarun.cell_ids==myCells(end));
    [X, Y] = drawEllipse([ctr(datInd,:) rad(datInd,:) fit_angle(datInd)]);
    [X, Y] = tformfwd(coord_tform, X, Y);
    hold on
    plot(X,Y,'b')
    
    deleted_cells_list=[deleted_cells_list myCells(end)];
    deleted_cells_list = unique(deleted_cells_list);
    deleted_cells_list = setdiff(deleted_cells_list, myCells(1:end-1));
    
    % update this cell info
    
    deletedCellUI=uicontrol('style','text', 'Units', 'Normalized','position',[0.4 0.7 0.08 0.2],...
        'string',{'',['Deleted cells ',int2str(deleted_cells_list)]},...
        'BackgroundColor',[.3 .5 .3],'ForegroundColor','k','fontsize',16);
    
    % delete from arrays
    myCells(end)=[];
    cones(end)=[];
    
    % update all cones plot
    show_all_cones;
    
    % update cone info
    update_cone_info;
    
elseif ~isempty(deleted_cells_list) % new cell added; checks if it was deleted
    
    deleted_cells_list = setdiff(deleted_cells_list, myCells);
    if ~isempty(deleted_cells_list) % if it was not the only cell
        deletedCellUI=uicontrol('style','text', 'Units', 'Normalized','position',[0.4 0.7 0.08 0.2],...
            'string',{'',['Deleted cells ',int2str(deleted_cells_list)]},...
            'BackgroundColor',[.3 .5 .3],'ForegroundColor','k','fontsize',16);
    else
        tmp=findobj('position',[0.4 0.7 0.08 0.2]);
        delete(tmp)
    end
end
        
