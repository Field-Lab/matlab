function add_cone(cone_list, contour_list, sta, sta_plot)

subplot(sta_plot);

while 1
    
    xlims = get(gca,'Xlim');
    ylims = get(gca,'Ylim');

    [x,y]=ginput(1);
    x = round(x);
    y = round(y); 
    if x<xlims(1) || x>xlims(2) || y<ylims(1) || y>ylims(2)
        break
    end
  
    tmp = sta(y-1:y+1,x-1:x+1);
    [~,ic]=sort(tmp(:), 'descend');
    [r, c]=ind2sub([3 3],ic(1:4));
    r = y + r -2;
    c = x + c -2;
    conts = plot_contour(r,c);
    plot(conts(:,1), conts(:,2),'r', 'linewidth',3)
    
    cone_list(end+1,1:2) = [x y];
    cone_list(end,3) = 1;    
    contour_list{end+1} = [r,c];
    
    current_active_cones = [current_active_cones size(cone_list,1)];
    
end
