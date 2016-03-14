function select_cells(all_sta, pols, rad, ind)

global cones cone_regions
persistent cone_peaks new_cone_regs


if ~exist('hmain', 'var')
    hmain=figure;
    set(hmain,'Name','Main Window','Position',[-1843  80 1810 1025],...
        'ToolBar','Figure', 'KeyPressfcn', @figure1_WindowKeyPressFcn, ...
        'WindowButtonDownFcn',@wbdcb)
end
hPlot=subplot('position',[0.25 0.05 0.72 0.9], 'SortMethod','childorder');

ind = ind-1;
radius = rad;
threshold = 3.5;
bord = 50;
tmp = size(all_sta);

allPlot=subplot('position',[0.02 0.6 0.21 0.3], 'SortMethod','childorder');
% set(allPlot, 'dataaspectRatio', [1 1 1])
axis([0, tmp(2), 0, tmp(1)])
axis ij
hold on

conePlot=subplot('position',[0.02 0.2 0.21 0.3], 'SortMethod','childorder');


ncells = tmp(end);
if size(tmp,2)==4
    rgb = 1;
else
    rgb = 0;
end
sta2show = [];
sta = [];
new_cone_color = [1 0 0];

[contour, map_speck] = contruct_cone_region(radius);


    function figure1_WindowKeyPressFcn(hmain, eventdata, handles)
        
        switch eventdata.Key
            case 'space'
                
                if ind<ncells
                    ind = ind+1;
                    if ~isempty(cones)
                        cones = [cones; cone_peaks];
                        cone_regions = [cone_regions new_cone_regs];
                        if cones(1)==0
                            cones(1,:) = [];
                        end
                        subplot(allPlot)
                        if ~isempty(new_cone_regs)
                            for i = 1:length(new_cone_regs)
                                plot(new_cone_regs{i}(:,1), new_cone_regs{i}(:,2));
                            end
                        else
                            for i = 1:length(cone_regions)
                                plot(cone_regions{i}(:,1), cone_regions{i}(:,2));
                            end
                        end
                            
                        subplot(hPlot)
                    else
                        cones = [0 0];
                        cone_regions = [];
                    end
                    
                    if rgb
                        sta=all_sta(:,:,:,ind);
                        sta2show = 1/max(abs(sta(:)));
                        sta2show = sta*sta2show/2+0.5;
                    else
                        sta=all_sta(:,:,ind);
                        sta2show = sta;
                        

                        
                    end
                    
                    
                    subplot(hPlot);
                    hold off
                    colormap gray
                    imagesc(sta2show);
                    set(gca, 'DataAspectRatio', [1 1 1])
                    hold on
                    if size(cones,1)== size(cone_regions,2)
                        title(['Cell ', int2str(ind), ' out of ', int2str(ncells), '; ', ...
                            int2str(size(cones, 1)), ' cones']);
                    elseif ind==1
                        title(['Cell ', int2str(ind), ' out of ', int2str(ncells)]);
                    else
                        subplot(allPlot)
                        text(800, 300, ['Cell ', int2str(ind), ' out of ', int2str(ncells), '; CHECK SIZES'], ...
                            'color', 'r', 'fontsize', 25);
                    end
                    
                    statmp = pols(ind)*sta;
                    
                    [cons, cone_peaks] = find_cones_auto(statmp, threshold);
                    cone_peaks = round(cone_peaks);
                    tmp = find(cone_peaks(:,2)>size(sta,1)-radius*2 | cone_peaks(:,2)<radius*2);
                    tmp = [tmp; find(cone_peaks(:,1)>size(sta,2)-radius*2 | cone_peaks(:,1)<radius*2)];
                    cone_peaks(tmp,:) = [];
                    
                    tmp = triu(pdist2(cone_peaks, cone_peaks));
                    tmp(tmp==0) = 100;
                    
                    [tmp, ~] = find(tmp<radius*2);
                    if ~isempty(tmp)
                        cone_peaks(tmp,:) = [];
                    end
                                        
                    tmp = pdist2(cones, cone_peaks);
                    tmp = min(tmp);
                    tmp = find(tmp<radius);
                    if ~isempty(tmp)
                        cone_peaks(tmp,:) = [];
                    end

                    if ~isempty(cone_peaks(:))
                        
                      
                        
                        close_cones = cones(:,1)>min(cone_peaks(:,1))-bord*2 & ...
                            cones(:,1)<max(cone_peaks(:,1))+bord*2 & ...
                            cones(:,2)>min(cone_peaks(:,2))-bord*2 & ...
                            cones(:,2)<max(cone_peaks(:,2))+bord*2;
                        
                        for i=find(close_cones)'
                            plot(cone_regions{i}(:,1), cone_regions{i}(:,2), 'color', [1 1 0]);
                        end
                        
                        
                        new_cone_regs = cell(1,size(cone_peaks,1));
                        coll_sta = zeros(size(all_sta,1),size(all_sta,2) );
                        for i=1:size(cone_peaks,1)
                            raw_row = map_speck(:,1)+cone_peaks(i,1);
                            raw_col = map_speck(:,2)+cone_peaks(i,2);
                            fin = [];row_ind = 1;
                            for row_adjust=-1:1
                                col_ind = 1;
                                for col_adjust = -1:1
                                    radj = raw_row + row_adjust;
                                    cadj = raw_col + col_adjust;
                                    linearInd = sub2ind([size(all_sta,1),size(all_sta,2)], cadj, radj);
                                    map = zeros(size(all_sta,1),size(all_sta,2));
                                    map(linearInd) = 1;
                                    if rgb
                                        fin(row_ind,col_ind) = sum(sum(sum(statmp.*repmat(map,1,1,3))));
                                    else
                                        fin(row_ind,col_ind) = sum(sum(statmp.*map));
                                    end
                                    col_ind = col_ind+1;
                                end
                                row_ind = row_ind+1;
                            end
                            [row_ind,col_ind] = find(fin ==max(fin(:)),1);
                            
                            tmp = -1:1;
                            cone_peaks(i,1) = cone_peaks(i,1) + tmp(row_ind);
                            cone_peaks(i,2) = cone_peaks(i,2) + tmp(col_ind);
                            new_cone_regs{i} = [contour(:,1)+cone_peaks(i,1) contour(:,2)+cone_peaks(i,2)];
                            
                            if pols(ind)>0;
                                col = 'g';
                            else
                                col = 'b';
                            end
                            subplot(hPlot)
                            plot(new_cone_regs{i}(:,1),new_cone_regs{i}(:,2), col);
                             
                        end
                        
                        axis([min(cone_peaks(:,1))-bord ...
                            max(cone_peaks(:,1))+bord ...
                            min(cone_peaks(:,2))-bord ...
                            max(cone_peaks(:,2))+bord]);
                        
                        subplot(conePlot)
                        imagesc(sta2show)
 
                    end
                else
                    cones = [cones; cone_peaks];
                    cone_regions = [cone_regions new_cone_regs];
                    delete(hmain)
                end
                
                
            case 'x' % delete entire cell
                cone_peaks = [];
                new_cone_regs =[];
                delete(hPlot.Children(1:end-1))
                
            case 'leftarrow'
                radius = radius-0.5;
                radius
                [contour, map_speck] = contruct_cone_region(radius);
                
            case 'rightarrow'
                radius = radius+0.5;
                radius
                [contour, map_speck] = contruct_cone_region(radius);
                
            case 'a' % find more cones: decrease threshold
                
                threshold = threshold-0.2;
                threshold
                
                subplot(hPlot);
                hold off
                colormap gray
                imagesc(sta2show);
                set(gca, 'DataAspectRatio', [1 1 1])
                hold on
                
                
                [cons, cone_peaks] = find_cones_auto(pols(ind)*sta, threshold);
                cone_peaks = round(cone_peaks);
                new_cone_regs = cell(1,size(cone_peaks,1));
                for i=1:size(cone_peaks,1)
                    new_cone_regs{i} = [contour(:,1)+cone_peaks(i,1) contour(:,2)+cone_peaks(i,2)];
                    if pols(ind)>0;
                        col = 'g';
                    else
                        col = 'b';
                    end
                    plot(new_cone_regs{i}(:,1),new_cone_regs{i}(:,2), col);
                end
                
                axis([min(cone_peaks(:,1))-bord ...
                    max(cone_peaks(:,1))+bord ...
                    min(cone_peaks(:,2))-bord ...
                    max(cone_peaks(:,2))+bord]);
                
                
            case 'z' % find fewer cones: increase threshold
                
                threshold = threshold+0.2;
                threshold
                
                subplot(hPlot);
                hold off
                colormap gray
                imagesc(sta2show);
                set(gca, 'DataAspectRatio', [1 1 1])
                hold on
                
                [cons, cone_peaks] = find_cones_auto(pols(ind)*sta, threshold);
                cone_peaks = round(cone_peaks);
                
                new_cone_regs = cell(1,size(cone_peaks,1));
                for i=1:size(cone_peaks,1)
                    new_cone_regs{i} = [contour(:,1)+cone_peaks(i,1) contour(:,2)+cone_peaks(i,2)];
                    if pols(ind)>0;
                        col = 'g';
                    else
                        col = 'b';
                    end
                    plot(new_cone_regs{i}(:,1),new_cone_regs{i}(:,2), col);
                end
                
                axis([min(cone_peaks(:,1))-bord ...
                    max(cone_peaks(:,1))+bord ...
                    min(cone_peaks(:,2))-bord ...
                    max(cone_peaks(:,2))+bord]);
                
        end
        
    end



    function wbdcb(src,callbackdata)
        
        seltype = src.SelectionType;
        
        % new location
        cp = hPlot.CurrentPoint;
        xinit = round(cp(1,1));
        yinit = round(cp(1,2));
        
        if strcmp(seltype,'normal') % left button ADD cone
            cone_peaks(end+1,:) = [xinit yinit];
            new_cone_regs{end+1} = [contour(:,1)+cone_peaks(end,1) contour(:,2)+cone_peaks(end,2)];
            plot(new_cone_regs{end}(:,1), new_cone_regs{end}(:,2), 'color', new_cone_color);
            
            axis([min(cone_peaks(:,1))-bord ...
                max(cone_peaks(:,1))+bord ...
                min(cone_peaks(:,2))-bord ...
                max(cone_peaks(:,2))+bord]);
            
            
        elseif strcmp(seltype,'alt') % right button DELETE cone
            
            [~,t] = min(pdist2([xinit yinit], cone_peaks));
            a = findobj('XData', new_cone_regs{t}(:,1), 'YData', new_cone_regs{t}(:,2));
            for i=1:length(a)
                delete(a(i));
            end
            new_cone_regs(t) = [];
            cone_peaks(t,:) = [];
            
            axis([min(cone_peaks(:,1))-bord ...
                max(cone_peaks(:,1))+bord ...
                min(cone_peaks(:,2))-bord ...
                max(cone_peaks(:,2))+bord]);
            
        elseif strcmp(seltype,'extend') % double click any button DELETE OLD cone
            
            [~,t] = min(pdist2([xinit yinit], cones));
            a = findobj('XData', cone_regions{t}(:,1), 'YData', cone_regions{t}(:,2));
            for i=1:length(a)
                delete(a(i));
            end
            cone_regions(t) = [];
            cones(t,:) = [];
            
            
        end
    end
end
