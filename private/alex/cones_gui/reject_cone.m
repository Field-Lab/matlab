function reject_cone(sta_plot)

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

    [~, ind]=min(pdist2([x y], cone_list(current_active_cones, 1:2)));

    cone_list(end,3) = 1;    
    contour_list{end+1} = [r,c];
    
end
