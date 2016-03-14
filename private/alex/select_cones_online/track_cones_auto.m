function [new_cones, new_regions] = track_cones_auto(all_sta, marks, pols, cones, radius)


[contour, map_speck] = contruct_cone_region(radius);
field_size = [size(all_sta,1), size(all_sta,2)];


t = ceil(length(cones)/20);
disp('********************')

for i=1:length(cones)
    if mod(i,t)==0
        fprintf('*')
    end
    
    x = round(cones(i,1));
    y = round(cones(i,2));
    
    a = [];
    for j=-1:1
        for k=-1:1
            a = [a find(marks(x-j,y-k,:))'];
        end
    end
    a = unique(a);
    
    if ~isempty(a) && y<field_size(1) && x<field_size(2)
        statmp = 0;
        for j=1:length(a)
            sta=pols(a(j))*all_sta(:,:,a(j));
            statmp = statmp+sta;
        end
       
        
        
        raw_row = map_speck(:,1)+x;
        raw_col = map_speck(:,2)+y;
        fin = [];row_ind = 1;
        for row_adjust=-1:1
            col_ind = 1;
            for col_adjust = -1:1
                radj = raw_row + row_adjust;
                cadj = raw_col + col_adjust;
                linearInd = sub2ind(field_size, cadj, radj);
                map = zeros(field_size);
                map(linearInd) = 1;
%                 if rgb
%                     fin(row_ind,col_ind) = sum(sum(sum(statmp.*repmat(map,1,1,3))));
%                 else
                    fin(row_ind,col_ind) = sum(sum(statmp.*map));
%                 end
                col_ind = col_ind+1;
            end
            row_ind = row_ind+1;
        end
        [row_ind,col_ind] = find(fin ==max(fin(:)),1);
        
        tmp = -1:1;
        
                
        new_cones(i,:) = [x + tmp(row_ind),  y + tmp(col_ind)];
        new_regions{i} = [contour(:,1)+new_cones(i,1) contour(:,2)+new_cones(i,2)];
                
        
              
        
        %         % get new cone position
        %         bord = 2;
        %         tmp = sta2show(y-bord:y+bord, x-bord:x+bord,:);
        %         tmp = imresize(tmp, 3);
        %         tmp = tmp(4:12,4:12,:);
        %
        %         xcoord = -bord*2:bord*2;
        %         xcoord = repmat(xcoord,bord*2*2+1,1);
        %         xcoord = repmat(xcoord,1,1,3);
        %         xcoord = tmp.*xcoord;
        %
        %         ycoord = [-bord*2:bord*2]';
        %         ycoord = repmat(ycoord,1, bord*2*2+1);
        %         ycoord = repmat(ycoord,1,1,3);
        %         ycoord = tmp.*ycoord;
        %
        %         new_x = sum(xcoord(:))/sum(tmp(:))/3;
        %         new_y = sum(ycoord(:))/sum(tmp(:))/3;
        %
        %         th = 0:pi/250:2*pi;
        %         x = cones(i,1)+new_x;
        %         y = cones(i,2)+new_y;
        %
        %         xunit = round(((radius * cos(th) + x)*4))/4;
        %         yunit = round(((radius * sin(th) + y)*4))/4;
        %
        %         a = [xunit; yunit]';
        %         m = [];
        %         for j=1:499
        %             if xunit(j)==xunit(j+1) && xunit(j)==xunit(j+2) && yunit(j)==yunit(j+1) && yunit(j+1)==yunit(j+2)
        %                 m = [m j];
        %             end
        %         end
        %         a(m,:) = [];

        
    end
end
