function [new_cones, new_regions] = track_cones_auto(all_sta, marks, pols, cones, radius)

x = 0;
y = 0;
th = 0:pi/250:2*pi;
xunit = round(((radius * cos(th) + x)*2))/2;
yunit = round(((radius * sin(th) + y)*2))/2;

a = [xunit; yunit]';
m = [];
for j=1:499
    if xunit(j)==xunit(j+1) && xunit(j)==xunit(j+2) && yunit(j)==yunit(j+1) && yunit(j+1)==yunit(j+2)
        m = [m j];
    end
end
a(m,:) = [];

[X, Y] = meshgrid(linspace(x-10,x+10, 21), linspace(y-10,y+10, 21));
in = inpolygon(X, Y, a(:,1), a(:,2));
[rows, cols] = find(in);
rows = rows-11;
cols = cols-11;

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
    
    if ~isempty(a) && y<297 && x<397
        tmp_sta = 0;
        for j=1:length(a)
            sta=pols(a(j))*all_sta(:,:,:,a(j));
            tmp_sta = tmp_sta+sta;
        end
        sta2show = tmp_sta/max(tmp_sta(:));
        sta2show(sta2show<0) = 0;
        sta2show = imresize(sta2show, 2, 'method', 'nearest');
        
        fin = [];p2 = 1;
        for k=-2:1:2
            p1 = 1;
            for kk = -2:1:2
                x = round(cones(i,1)*2+k);
                y = round(cones(i,2)*2+kk);
                
                trows = rows+x;
                tcols = cols+y;
                linearInd = sub2ind([600,800], tcols, trows);
                map = zeros(600,800);
                map(linearInd) =1;
                
                fin(p1,p2) = sum(sum(sum(sta2show.*repmat(map,1,1,3))));
                p1 = p1+1;
            end
            p2 = p2+1;
        end
        [kk,k] = find(fin ==max(fin(:)),1);
        
        tmp = -2:1:2;
        x = round(cones(i,1)*2 + tmp(k));
        y = round(cones(i,2)*2 + tmp(kk));
        
        trows = rows+x;
        tcols = cols+y;
        linearInd = sub2ind([600,800], tcols,trows);
        map = uint8(zeros(600,800));
        map(linearInd) =1;
        
        dd = imresize(map,5,'method', 'nearest');
        [r, c] = find(dd,1);
        contour = bwtraceboundary(dd,[r c],'W',8,Inf,'counterclockwise');
        contour= round(contour/5)+0.5;
        
        tmp = [];
        for j=1:length(contour)-2
            if contour(j,1)==contour(j+1,1)  && contour(j,2)==contour(j+1,2)
                tmp = [tmp j];
            end
        end
        contour(tmp,:) = [];
%         
%         figure;
%         imagesc(sta2show)
%         hold on
%         plot(contour(:,1), contour(:,2), 'r')
        
         
        new_regions{i} = contour;
        new_cones(i,:) = [x y];
        
        
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
