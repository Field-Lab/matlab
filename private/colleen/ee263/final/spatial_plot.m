function [] = spatial_plot(x , y , z , nlevels , cmap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spatial_plot(x , y , z , nlevels , cmap)
% makes a scatter plot of the points (x(i),y(i)), where
% each point is color-coded to indicate the intensity z(i);
% the color codes contain nlevels, and the optional cmap
% argument can be used to specify a color map; if cmap is
% not specified, the cool colormap is used
% -- x : a vector containing the x-coordinates
% -- y : a vector containing the y-coordinates
% -- z : a vector containing the intensity
% -- nlevels : the number of intensity levels to use
% -- cmap : the color map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check the number of arguments, and use the default
    % color map if no color map was provided
    if nargin < 4
        error('not enough arguments provided to spatial_plot');
    end
    if nargin == 4
        cmap = cool(nlevels);
    end

    % compute the levels
    zmin = min(z);
    zmax = max(z);
    delta = (zmax - zmin + 1e-6) / nlevels;
    
    % make a color-coded scatter plot of the data points
    figure();
    hold on;
    for i = 1:length(x)
        k = 1 + floor((z(i) - zmin) / delta);
        plot(x(i) , y(i) , 'o' , ...
             'Color' , cmap(k,:) , ...
             'MarkerFaceColor' , cmap(k,:));
    end
    hold off;
    box on;
    
    % set the limits of the axes
    xmin = min(x);
    xmax = max(x);
    xlim([xmin xmax] + 0.1 * (xmax - xmin) * [-1 +1]);

    ymin = min(y);
    ymax = max(y);
    ylim([ymin ymax] + 0.1 * (ymax - ymin) * [-1 +1]);
end