function show_next_sta(sta, current_cell, cone_list, sta_plot)

subplot(sta_plot);
colormap gray
hold off
imagesc(sta)
hold on

plot(cone_peaks(:,2), cone_peaks(:,1), 'xr')

% weights
t = [cone_peaks(:,2), cone_peaks(:,1)];
for j=1:length(cone_peaks)
    m(j) = sta(t(j,2),t(j,1));
    if m(j)<0; m(j) = 0; end
    plot(t(j,1), t(j,2), '+g', 'markersize', m(j)*100+0.1)
end
k = m/max(m);
kk = find(k>0.3);
col = 'bykrgbykrgbykrgbykrg';
for j=1:length(kk)
    tmp = sta(t(kk(j),2)-1:t(kk(j),2)+1,t(kk(j),1)-1:t(kk(j),1)+1);
    [~,ic]=sort(tmp(:), 'descend');
    [r, c]=ind2sub([3 3],ic(1:4));
    r = t(kk(j),2) + r -2;
    c = t(kk(j),1) + c -2;
    conts = plot_contour(r,c);
    plot(conts(:,1), conts(:,2),col(j), 'linewidth',3)
end

current_active_cones
    

