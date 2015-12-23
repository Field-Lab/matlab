% for each match, plot the axon, and a summary of the EI which might point to the soma


switch 1
    case 1
        % Show the image
        im = fstitch_arr;
        imshow(im, 'XData', fstitch_xd, 'YData', fstitch_yd);
        axis xy;
        hold on;
        axon_color = 'w';
    case 2
        im = a2tf;
        imshow(im, 'XData', a2xd, 'YData', a2yd);
        axis xy;
        hold on;
        axon_color = 'w';
    case 3
        hold on;
        axon_color = 'g';
end

% for each match, plot the axon and EI estimate
nlp = datarun.ei.nlPoints;
pos = datarun.ei.position;

for mm = 1:size(match_list,1)
    axon_id = match_list(mm,1);
    cell_id = match_list(mm,2);
    
    % plot the axon
    points = traced_cell_points(axons{axon_id}(1:2,:),axons{axon_id}(2:end,:));
    plot(points(:,1), points(:,2), axon_color);

%     text(axons{axon_id}(1,1),axons{axon_id}(1,2),sprintf('axon %d, EI %d',axon_id,cell_id),'Color','g','FontSize',15,...
%                 'HorizontalAlignment','Center','VerticalAlignment','Bottom')
     
     % get the EI com
     com = ei_com_(get_ei(datarun,cell_id),pos,nlp,'frames',[-8 0 24 -8:4],'roi',{'peak',60*3});
     com_line(1,:) = com{1}(1,:);  % -4, negative
     com_line(2,:) = com{1}(2,:);  %  0, negative
     com_line(3,:) = com{1}(3,:);  % 20, negative
     com_line(4,:) = com{2}(3,:);  % 20, positive
     com_line(5:17,:) = com{1}(4:16,:);  % 20, negative
     % plot it
     plot(com_line([1 2 4],1),com_line([1 2 4],2),'k')
     %plot(com_line([3 4],1),com_line([3 4],2),'m')
     plot(com_line(1,1),com_line(1,2),'.','color',[.5 .5 1],'MarkerSize',40)
     plot(com_line(2,1),com_line(2,2),'.','color',[.5 .5 1],'MarkerSize',20)
     %plot(com_line(3,1),com_line(3,2),'.','color',[.5 .5 1],'MarkerSize',10)
     plot(com_line(4,1),com_line(4,2),'.','color',[1 .5 .5],'MarkerSize',10)
     plot(com_line(5:17,1),com_line(5:17,2),'b.-.')
     
     if axon_id == 65; disp(com_line); end
     
     % plot line connecting them
%     plot([com_line(4,1) axons{axon_id}(1,1)],[com_line(4,2) axons{axon_id}(1,2)],'color',[.8 .8 .8])
     
     % add line to largest electrode
%     ei = get_ei(datarun,cell_id);
%     [junk,peak] = max(max(ei,[],2));
%     plot([pos(peak,1) axons{axon_id}(1,1)],[pos(peak,2) axons{axon_id}(1,2)],'color','m')
end

% add electrodes
%plot(pos(:,1),pos(:,2),'o','color',[.5 .5 0])

% set xlim, ylim
axis equal
%xlim(datarun.ei.array_bounds_x)
%ylim(datarun.ei.array_bounds_y)
