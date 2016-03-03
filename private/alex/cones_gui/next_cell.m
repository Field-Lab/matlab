function next_cell(frame)

global datarun current_cell


while 1
    subplot(hplot);
    plot_position=get(hplot,'position');
    [x,y]=ginput(1);
    clicked_plot_position=get(gca,'position');
    
    dd=nnz(plot_position-clicked_plot_position);
    
    if dd || x<1 || x>datarun.stimulus.field_width || y<1 || y>datarun.stimulus.field_height
        return;
    end
        
    tmp=pdist2([x,y],ctr);
    [~,ind]=min(tmp(cell_indices));    
    [X, Y] = drawEllipse([ctr(cell_indices(ind),:) rad(cell_indices(ind),:) fit_angle(cell_indices(ind))]);
    [X, Y] = tformfwd(coord_tform, X, Y);
    IN = inpolygon(x,y,X,Y);
    if IN % clicked inside
        
        datInd=datarun.cell_ids(cell_indices(ind));
        
        % check if this cell has been selected already
        if ~isempty(find(myCells==datInd, 1))
            
            % make this cell the last in the list (for convenience in add manually)
            tmp=cones{myCells==datInd};
            cones(myCells==datInd)=[];
            cones{end+1}=tmp;            
            myCells(myCells==datInd)=[];
            myCells=[myCells datInd];
            
            show_sta(cell_indices(ind),'index_type','datarun_id','frame', frame);
            show_cones(1,cell_indices(ind),'index_type','datarun_id');
            return;
        else
            
            % color rf fit on main plot
            plot(X,Y,'r')
            myCells=[myCells datInd];
            show_sta(cell_indices(ind),'index_type','datarun_id','frame', frame);
            handle_cones(0,cell_indices(ind),'index_type','datarun_id', 'frame', frame, 'nnd_scale', nnd_scale)
            
            % add more cones
            flag=1;
            while flag
                choice = questdlg('Add a cone?', 'Add a cone', ...
                    'Yes please!','Kill it with fire','No','Yes please!');
                % Handle response
                switch choice
                    case 'Yes please!'
                        handle_cones(0,cell_indices(ind),'index_type','datarun_id')
                        
                    case 'Kill it with fire'
                        handle_cones(3,cell_indices(ind),'index_type','datarun_id')
                        flag=0;
                        
                    case 'No'
                        flag=0;
                end
            end
            
        end
        
    end
end
