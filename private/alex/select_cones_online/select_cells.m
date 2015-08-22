function select_cells(hplot, nnd_scale, frame)

global myCells datarun cones cell_indices coord_tform ctr rad fit_angle


while 1
    subplot(hplot);
    plot_position=get(hplot,'position');
    [x,y]=ginput(1);
    clicked_plot_position=get(gca,'position');
    
    dd=nnz(plot_position-clicked_plot_position);
    
    if dd || x<1 || x>datarun.stimulus.field_width || y<1 || y>datarun.stimulus.field_height
        return;
    else
        
    tmp=pdist2([x,y],ctr);
    [~,ind]=min(tmp(cell_indices));    
    [X, Y] = drawEllipse([ctr(cell_indices(ind),:) rad(cell_indices(ind),:) fit_angle(cell_indices(ind))]);
    [X, Y] = tformfwd(coord_tform, X, Y);
%     IN = inpolygon(x,y,X,Y);
%     if IN % clicked inside
        
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
%         
%         
%         myInd=datarun.cell_ids(cell_indices(ind));
%         datarun = load_sta(datarun, struct('load_sta',myInd));
%         
%         
%         % plot sta
%         if ishandle(hStaPlot)
%             delete(hStaPlot)
%         end
%    
%         hStaPlot=subplot('position',[0.4 0.62 0.3 0.3]);       
%         set(gca,'DataAspectRatio',[1 1 1])
%         axis ij
%         hold on
%         
%         centre=round(ctr(cell_indices(ind),:));
%         radius=ceil(max(rad(cell_indices(ind),:)))*scale_factor;
%         
%         sta=squeeze(datarun.stas.stas{cell_indices(ind)});  
%         sta=sta(:,:,5);
%         colormap gray
%         imagesc(sta)
%         
%         plot(X,Y,'g')
%         
%         minx=max(centre(1)-radius,1);
%         maxx=min(centre(1)+radius,datarun.stimulus.field_width);
%         miny=max(centre(2)-radius,1);
%         maxy=min(centre(2)+radius, datarun.stimulus.field_height);
%         
%        
%         
%         % find cones
%         tmp_sta=sta;
%         if abs(min(tmp_sta(:)))>max(tmp_sta(:))  % OFF cell, invert polarity
%             tmp_sta=-tmp_sta;
%         end
%         
%         w_center=[];
%         max_dist=max(rad(cell_indices(ind),:))*3;
%         cones{end+1}=[];
%         cnt=0;
%         flag=1;
%         while flag
%             myMax=max(tmp_sta(:));
%             myCoord=find(tmp_sta==myMax,1);
%             [x, y]=ind2sub(size(tmp_sta),myCoord);
%             % check if border position... if, then correct - don't trust!
%             x=min(max(x,3),size(sta,1)-3);
%             y=min(max(y,3),size(sta,2)-3);
%             
%             if cnt>4 % start checking from 6th cone, 5 are obligatory
%                 % center of mass
%                 x_center=cones{end}(:,1);
%                 y_center=cones{end}(:,2);
%                 x_com=double(sum(w_center.*x_center)/sum(w_center));
%                 y_com=double(sum(w_center.*y_center)/sum(w_center));
%                 if pdist2([y x], [x_com y_com])>3*max_dist
%                     adjust_axis(cones, minx, maxx, miny, maxy)
%                     flag=0;
%                 end
%                 
%                 % mean nnd
%                 tmp=squareform(pdist([cones{end}(:,1),cones{end}(:,2)]));
%                 tmp(tmp==0)=10000;
%                 tmp=min(tmp);
%                 if min(pdist2([y x], [cones{end}(:,1),cones{end}(:,2)]))>(mean(tmp)+scale_distance*std(tmp))
%                     adjust_axis(cones, minx, maxx, miny, maxy)
%                     flag=0;
%                 end
%             end 
%             if flag
%                 cnt=cnt+1;
%                 w_center=[w_center; tmp_sta(x,y)];
%                 tmp_sta(x-2:x+2,y-2:y+2,:)=0;
%                 plot(y,x,'x','color','r','markersize',20)
%                 cones{end}=[cones{end}; [y, x]];
%                 
%             end
%         end
%         if length(cones)>1
%             for k=1:length(cones)
%                 plot(cones{k}(1,:),cones{k}(2,:),'y+')
%             end
%         end
%         
%         % add cones one by one
%         flag=1;
%         while flag
%             choice = questdlg('Add a cone?', 'Add a cone', ...
%                 'Yes please!','Kill it with fire','No','Yes please!');
%             % Handle response
%             switch choice
%                 case 'Yes please!'
%                     myMax=max(tmp_sta(:));
%                     myCoord=find(tmp_sta==myMax,1);
%                     [x, y]=ind2sub(size(tmp_sta),myCoord);
%                     
%                     % check if border position... if, then correct - don't trust!
%                     x=min(max(x,3),size(sta,1)-3);
%                     y=min(max(y,3),size(sta,2)-3);
%                     
%                     tmp_sta(x-2:x+2,y-2:y+2,:)=0;
%                     plot(y,x,'x','color','r','markersize',20)
%                     cones{end}=[cones{end}; [y, x]];     
%                     adjust_axis(cones, minx, maxx, miny, maxy)
% 
%                     cnt=cnt+1;
%                 case 'Kill it with fire'
%                     tmp_sta(x-2:x+2,y-2:y+2,:)=sta(x-2:x+2,y-2:y+2,:);
%                     
%                     children = get(gca, 'children');
%                     delete(children(1));
%                     
%                     cones{end}=cones{end}(1:end-1,:);
%                     adjust_axis(cones, minx, maxx, miny, maxy)
% 
%                     cnt=cnt-1;
%                     flag=0;
%                 case 'No'
%                     flag=0;
%             end
%         end
%         
%     end
% end
