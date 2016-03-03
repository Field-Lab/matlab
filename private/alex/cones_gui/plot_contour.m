function conts = plot_contour(rows,cols)

t = zeros(max(rows)+1, max(cols)+1);

for i=1:length(rows)
    t(rows(i), cols(i)) = 1;
end

dd = imresize(t,5,'method', 'nearest');
[r, c] = find(dd,1);
contour = bwtraceboundary(dd,[r c],'W',8,Inf,'counterclockwise');
contour= round(contour/5)+0.5;

conts = [contour(:,2),contour(:,1)];
