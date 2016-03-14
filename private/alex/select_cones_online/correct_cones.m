function correct_cones(weights, all_sta, pols, cones, cone_regions, cone_lists, radius, weight_threshold)

global new_cones new_regions
new_cones = [];
new_regions = [];
i = 0;
hmain = figure;
set(hmain, 'position', [-1650  584  1130 500], 'KeyPressfcn', @figure1_WindowKeyPressFcn, ...
    'WindowButtonDownFcn',@wbdcb)
tPlot=subplot(1,2,1);
hPlot=subplot(1,2,2);
[contour, map_speck] = contruct_cone_region(radius);

colors = 'rkrkrkrkrkrkrrkrkrkrkrkrkrrkrkrkrkrkrkrrkrkrkrkrkrkrrkrkrkrkrkrkr';

    function figure1_WindowKeyPressFcn(hmain, eventdata, handles)
        
        switch eventdata.Key
            case 'space'
                if i<length(cone_lists)
                    i = i+1
                    
                    suspect_cones = cone_lists{i};
                    % find cells with significant weight on the cones
                    a = [];
                    for j=1:length(suspect_cones)
                        a = [a find(weights(suspect_cones(j), :)>weight_threshold)];
                    end
                    
                    a = unique(a);
                    
                    % max 3 cells... if more overlap, repeat 2 times - shouldn't be a big
                    % problem
                    comb = 0;
                    for j=1:length(a)
                        b = pols(a(j))*all_sta(:,:,a(j));
                        comb =  comb + b/max(b(:));
                    end
                    comb = comb/length(a);
                    
%                     for j=1:min(3, length(a))
%                         b = pols(a(j))*all_sta(:,:,a(j));
%                         comb(:,:,j) = b/max(b(:));
%                     end
                    
                    
                    subplot(tPlot)
                    hold off
                    imagesc(comb);
                    hold on
                    for j=1:length(suspect_cones)
                        plot(cone_regions{suspect_cones(j)}(:,1), cone_regions{suspect_cones(j)}(:,2), 'color', colors(j));
                    end
%                     plot(cone_regions{c(i)}(:,1), cone_regions{c(i)}(:,2));
                    axis([cones(suspect_cones(j),1)-10, cones(suspect_cones(j),1)+10,...
                        cones(suspect_cones(j),2)-10, cones(suspect_cones(j),2)+10])
                    set(gca, 'DataAspectRatio', [ 1 1 1])
                    title(int2str(length(a)))
                    
                    subplot(hPlot);
                    hold off
                    imagesc(comb);
                    hold on
                    axis([cones(suspect_cones(j),1)-10, cones(suspect_cones(j),1)+10,...
                        cones(suspect_cones(j),2)-10, cones(suspect_cones(j),2)+10])
                    set(gca, 'DataAspectRatio', [ 1 1 1])
                else
                    delete(hmain)
                end
        end
        
    end

    function wbdcb(src,callbackdata)
        
        seltype = src.SelectionType;
        
        % new location
        cp = hPlot.CurrentPoint;
        xinit = cp(1,1);
        yinit = cp(1,2);
        
        if strcmp(seltype,'normal') % left button ADD cone
            
            new_cones(end+1,:) = [xinit yinit];
            new_regions{end+1} = [contour(:,1)+new_cones(end,1) contour(:,2)+new_cones(end,2)];
            
            subplot(hPlot);
            plot(new_regions{end}(:,1), new_regions{end}(:,2), 'color', [1 1 1]);
            
        elseif strcmp(seltype,'alt') % right button DELETE cone
            
            [~,t] = min(pdist2([xinit yinit], new_cones));
            tmp = findobj('XData', new_regions{t}(:,1), 'YData', new_regions{t}(:,2));
            if ~isempty(tmp)
                delete(tmp(1));
                new_regions(t) = [];
                new_cones(t,:) = [];
            else
                disp('No new cones?')
                i
            end
        end
    end
end