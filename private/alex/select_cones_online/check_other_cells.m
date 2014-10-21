function check_other_cells(frame)

global datarun cones myCells ctr cell_indices

myCones=cell2mat(cones');
all_cones=zeros(5,5,length(cell2mat(cones')));
cnt=1;
for i=1:length(myCells)
    datInd=find(datarun.cell_ids==myCells(i));
    sta=squeeze(datarun.stas.stas{datInd});
    sta=sta(:,:,frame); % frame
    
    for j=1:length(cones{i})
        all_cones(:,:,cnt)=sta(cones{i}(j,2)-2:cones{i}(j,2)+2,cones{i}(j,1)-2:cones{i}(j,1)+2);
        cnt=cnt+1;
    end
    conesPerCell(i)=size(cones{i},1);
end

myArea=convhull(myCones(:,1),myCones(:,2));

figure
plot(myCones(:,1),myCones(:,2),'*r')
hold on
plot(myCones(myArea,1),myCones(myArea,2),'r')

myColors='brkcgm';
myMarkers='+xdv><';
cnt=1;
%find cells with RFs within considered area
for i=1:length(datarun.cell_ids)
    
    if isempty(find(cell_indices==i, 1)) && isempty(find(datarun.cell_types{1, end}.cell_ids==datarun.cell_ids(i), 1))
        IN = inpolygon(ctr(i,1),ctr(i,2),myCones(myArea,1),myCones(myArea,2));

        if IN
            % check and load sta if not loaded
            if isempty(datarun.stas.stas{i})
                datarun = load_sta(datarun, struct('load_sta',datarun.cell_ids(i)));
            end
            
            % prepare sta
            sta=squeeze(datarun.stas.stas{i});
            sta=sta(:,:,frame); % frame
            
            
            
            tmp_sta=sta;
            if abs(min(tmp_sta(:)))>max(tmp_sta(:))  % OFF cell, invert polarity
                tmp_sta=-tmp_sta;
            end
            
            weight_threshold=robust_mean(tmp_sta(:))+3*robust_std(tmp_sta(:));
            
            w_center=[];
            cones_tmp=[];
            
            keep_looking=true;
            while keep_looking
                myMax=max(tmp_sta(:));
                myCoord=find(tmp_sta==myMax,1);
                [row, col]=ind2sub(size(tmp_sta),myCoord);
                % check if border position... if, then correct - don't trust!
                row=min(max(row,3),size(sta,1)-3);
                col=min(max(col,3),size(sta,2)-3);
                
                if tmp_sta(row,col)>weight_threshold
                    
                    
                    if length(cones_tmp)>5 % start checking from 6th cone, 5 are obligatory
                        
                        % center of mass: check if new cone is not at most 2
                        % of max distance between center and other cones
                        x_center=cones_tmp(:,1);
                        y_center=cones_tmp(:,2);
                        x_com=double(sum(w_center.*x_center)/sum(w_center));
                        y_com=double(sum(w_center.*y_center)/sum(w_center));
                        
                        all_dist=pdist2([x_com y_com], [x_center,y_center]);
                        max_dist=max(all_dist);
                        
                        if pdist2([col row], [x_com y_com])>1.5*max_dist
                            keep_looking=false;
                        end
                        
                        % mean nnd
                        tmp=squareform(pdist([x_center,y_center]));
                        tmp(tmp==0)=10000;
                        tmp=min(tmp);
                        if min(pdist2([col row], [x_center,y_center]))>(mean(tmp)+2*std(tmp))
                            keep_looking=false;
                        end
                        
                        % relative weight
                    end
                else
                    keep_looking=false;
                end
               
                if keep_looking
                    w_center=[w_center; tmp_sta(row,col)];
                    tmp_sta(row-2:row+2,col-2:col+2,:)=0;
                    cones_tmp=[cones_tmp; [col, row]];
                end
            end
            if cnt>length(myColors)
                myC='k';
                myM='x';
            else
                myC=myColors(cnt);
                myM=myMarkers(cnt);
            end
            plot(cones_tmp(:,1),cones_tmp(:,2),[myM,myC]);
            cnt=cnt+1;
%             myColors
%             figure
%             colormap gray
%             imagesc(sta)
%             hold on
%             plot(cones_tmp(:,1),cones_tmp(:,2),'*r')
%             
            
            
            
            
            
            
        end
    end
end

axis ij
axis([0 200 0 200])
title([int2str(cnt-1), ' cells found!'])



