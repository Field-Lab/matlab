function laplacian_smoothing_plot(y1,text1,y2,text2,y3,text3)
%produces an imagesc and surface plot for each y value given,
%keeping the color and axis scale constant across all plots,
%and places the plots in a subplot window.

%each y should be a 400 element column vector, representing a 20x20 grid
%setting y3 = 0 will only produce plots for y1 and y2.

%each node is connected to the node
%up, down, left, and right of it, unless
%there is a black line drawn between the cells in the plot

if y3 == 0 %we only plot for y1 and y2
    n = 2;
    y = [y1 y2];
    text = {text1 text2};
else %plot for all three
    n = 3;
    y = [y1 y2 y3];
    text = {text1 text2 text3};
end

cmin = -25; %the expected smallest value in the data, for axis scaling and coloring
cmax = 15; %the expected largest value in the data

figure
for i = 1:n
    subplot(2,n,i)
    %reshape the vector y into a 20x20 matrix for plotting
    imagesc(reshape(y(:,i),20,20))
    caxis([cmin cmax]) %set color scale
    title(text{i})

    %draw the lines which mark that cells are not connected with an edge
    hold on
    plot(10.5*ones(13,1),7.5:.5:13.5,'k','LineWidth',3)
    plot(4.5:.5:16.5,7.5*ones(25,1),'k','LineWidth',3)
    plot(4.5:.5:16.5,13.5*ones(25,1),'k','LineWidth',3)
    hold off
    
    subplot(2,n,i+n)
    %reshape the vector y into a 20x20 matrix for plotting
    surf(reshape(y(:,i),20,20))
    caxis([cmin cmax]) %set color scale
    axis([0 20 0 20 cmin cmax]) %set z-axis scale
    %title(text{i})
end