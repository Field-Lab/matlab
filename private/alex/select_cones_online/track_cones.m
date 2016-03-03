function track_cones(all_sta, marks, pols, cones, cone_regions, radius, start_cone)

global new_cones new_regions
new_cones = [];
new_regions = [];
i = start_cone-1;
hmain = figure;
set(hmain, 'position', [-1650         157        1161         927], 'KeyPressfcn', @figure1_WindowKeyPressFcn, ...
    'WindowButtonDownFcn',@wbdcb)
oldPosPlot = subplot(2,2,1);
oldNegPlot = subplot(2,2,2);
fullViewPlot = subplot(2,2,3);
oldNegCone = [];
oldPosCone = [];
xinit = [];
yinit = [];
new_cone_region = [];

    function figure1_WindowKeyPressFcn(hmain, eventdata, handles)
        
        switch eventdata.Key
            case 'space'
                if i<length(cones)
                    if ~isempty(xinit)
                        new_cones(end+1,:) = [xinit yinit];
                        new_regions{end+1} = new_cone_region;
                    end
                    
                    flag = 1;
                    
                    while flag
                        i = i+1;
                        
                        x = round(cones(i,1));
                        y = round(cones(i,2));
                        
                        a = [];
                        for j=-1:1
                            for k=-1:1
                                a = [a find(marks(x-j,y-k,:))'];
                            end
                        end
                        a = unique(a);
                        if ~isempty(a)
                            tmp_sta = 0;
                            for j=1:length(a)
                                sta=pols(a(j))*all_sta(:,:,:,a(j));
                                tmp_sta = tmp_sta+sta;
                            end
                            sta2show = tmp_sta/max(tmp_sta(:));
                            sta2show(sta2show<0) = 0;
                            sta2show = imresize(sta2show, 2, 'method', 'nearest');
                            
                            subplot(oldPosPlot)
                            hold off
                            imagesc(sta2show)
                            common_axis = [x*2-15, x*2+15, y*2-15, y*2+15];
                            axis(common_axis)
                            hold on
%                             oldPosCone = plot(cone_regions{i}(:,1), cone_regions{i}(:,2),  'color', [1 1 1]*0.5);
                            title(['Cone ', int2str(i),', ', int2str(length(a)), ' cells found'])
                            
                            sta2showNeg = -sta2show+1;
                            subplot(oldNegPlot)
                            hold off
                            imagesc(sta2showNeg)
                            axis(common_axis)
                            hold on
%                             oldNegCone = plot(cone_regions{i}(:,1), cone_regions{i}(:,2),  'color', [0 0 0]+0.5);
                            title(['Cone ', int2str(i),', ', int2str(length(a)), ' cells found'])
                            
                            
                            % get new cone position
%                             bord = 2;
%                             tmp = sta2show(y-bord:y+bord, x-bord:x+bord,:);
%                             tmp = imresize(tmp, 3, 'method', 'nearest');
%                             tmp = tmp(4:12,4:12,:);
%                             tmp = tmp.^2;
%                
%                             xcoord = -bord*2:bord*2;
%                             xcoord = repmat(xcoord,bord*2*2+1,1);
%                             xcoord = repmat(xcoord,1,1,3);
%                             xcoord = tmp.*xcoord;
% 
%                             ycoord = [-bord*2:bord*2]';
%                             ycoord = repmat(ycoord,1, bord*2*2+1);
%                             ycoord = repmat(ycoord,1,1,3);
%                             ycoord = tmp.*ycoord;
                            
%                             figure
%                             imagesc(tmp)
%                             hold on
%                             plot(a(:,1)-cones(i,1)+5, a(:,2)-cones(i,2)+5,  'color', [1 0 0], 'linewidth', 2);

%                             new_x = (sum(xcoord(:))/sum(tmp(:)))/(9/5);
%                             new_y = (sum(ycoord(:))/sum(tmp(:)))/(9/5);                            
                            
                            th = 0:pi/250:2*pi;
 
                            fin = [];p2 = 1;
                            for k=-2:0.2:2
                                p1 = 1;
                                k
                                for kk = -2:0.2:2                                    
                                    x = cones(i,1)*2+k;
                                    y = cones(i,2)*2+kk;
                                    %
                                    xunit = round(((radius * cos(th) + x)*4))/4;
                                    yunit = round(((radius * sin(th) + y)*4))/4;
                                    
                                    a = [xunit; yunit]';
                                    m = [];
                                    for j=1:499
                                        if xunit(j)==xunit(j+1) && xunit(j)==xunit(j+2) && yunit(j)==yunit(j+1) && yunit(j+1)==yunit(j+2)
                                            m = [m j];
                                        end
                                    end
                                    a(m,:) = [];
                                                                        
                                    map = zeros(600,800);

                                    t = round([x y]);
                                    cr = a;
                                    [X, Y] = meshgrid(linspace(t(1)-10,t(1)+10, 21), linspace(t(2)-10,t(2)+10, 21));
                                    [in, on] = inpolygon(X, Y, cr(:,1), cr(:,2));
                                    in = find(in);
                                    for cnt=1:length(in)
                                        map(Y(in(cnt)), X(in(cnt))) = 1;
                                    end
                                    
%                                     figure
%                                     imagesc(tmp)
%                                     hold on
%                                     plot(a(:,1), a(:,2),  'color', [1 0 0], 'linewidth', 2);
                                    
                                    fin(p1,p2) = sum(sum(sum(sta2show.*repmat(map,1,1,3))));
                                    p1 = p1+1;
                                end
                                p2 = p2+1;
                            end
                            [kk,k] = find(fin ==max(fin(:)),1);
                            
                            tmp = -2:0.2:2;
                            x = cones(i,1)*2 + tmp(k);
                            y = cones(i,2)*2 + tmp(kk);
                            
                            %
                            xunit = round(((radius * cos(th) + x)*4))/4;
                            yunit = round(((radius * sin(th) + y)*4))/4;
                            
                            a = [xunit; yunit]';
                            m = [];
                            for j=1:499
                                if xunit(j)==xunit(j+1) && xunit(j)==xunit(j+2) && yunit(j)==yunit(j+1) && yunit(j+1)==yunit(j+2)
                                    m = [m j];
                                end
                            end
                            a(m,:) = [];
                            new_cone_region = a;

                            subplot(oldPosPlot)
                            tmpPosCone = plot(a(:,1), a(:,2),  'color', [1 0 0], 'linewidth', 2);
                            subplot(oldNegPlot)
                            tmpNegCone = plot(a(:,1), a(:,2),  'color', [1 0 0], 'linewidth', 2);
                            
                             
                            map = zeros(600,800);
                            t = round([x y]);
                            cr = a;
                            [X, Y] = meshgrid(linspace(t(1)-10,t(1)+10, 21), linspace(t(2)-10,t(2)+10, 21));
                            [in, on] = inpolygon(X, Y, cr(:,1), cr(:,2));
                            in = find(in);
                            for cnt=1:length(in)
                                map(Y(in(cnt)), X(in(cnt))) = 0.8;
                            end
                            map = repmat(map,1,1,3);
                            map(:,:,2:3) = 0;
                            map = map+sta2show;
                            subplot(fullViewPlot)
                            imagesc(map)  
                            axis(common_axis)
                            
                            flag = 0;
                        end
                    end
                    
                else
                    delete(hmain)
                end
        end
        
    end

    function wbdcb(src,callbackdata)
        
        seltype = src.SelectionType;
        
        % new location
        selectedPlot = gca;
        cp = selectedPlot.CurrentPoint;
        xinit = cp(1,1);
        yinit = cp(1,2);
        
        if strcmp(seltype,'normal') % left button CHANGE cone
            
            % first delete cone regions from plot
            if ishandle(oldPosCone)
                delete(oldPosCone);
            end
            if ishandle(oldNegCone)
                delete(oldNegCone);
            end
            
            % then add new cones
            new_cones(end+1,:) = [xinit yinit];
            th = 0:pi/250:2*pi;
            xunit = round(((radius * cos(th) + xinit)*4))/4;
            yunit = round(((radius * sin(th) + yinit)*4))/4;
            aa = [xunit; yunit]';
            m = [];
            for k=1:499
                if xunit(k)==xunit(k+1) && xunit(k)==xunit(k+2) && yunit(k)==yunit(k+1) && yunit(k+1)==yunit(k+2)
                    m = [m k];
                end
            end
            aa(m,:) = [];
            
            subplot(oldPosPlot);
            oldPosCone = plot(aa(:,1), aa(:,2), 'color', [1 1 1]);
            subplot(oldNegPlot);
            oldNegCone = plot(aa(:,1), aa(:,2), 'color', [0 0 0]);
            new_cone_region = aa;
            
        end
    end
end


