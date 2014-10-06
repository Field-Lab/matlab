function p = toScreen(p, figNum, electrode, clusters)
%TOSCREEN Plot.toScreen method
%   p = TOSCREEN(Plot p, figureNumber, electrodeNumber) 
%   generates a new plot window using figure number f and plots 
%   the output of the clustering algorithm on electrode e to the screen.
%
%   tamachado@salk.edu 1/23/08

% Get the data we need to plot
data = toStruct(p.data, electrode);

% Figure out how many clusters there are (in the MAP only, if multiple samples)
if nargin > 3
   nClusters = clusters.nClusters(1);
else
   nClusters = 1;
end


% Plot the pc projections with a density gradient
f.axes = generateAxes(figNum, nClusters);
f.image = densityGradient(data, f.axes.projection_axes);
p.figure{electrode} = f;


% If clusters are provided, plot them
if nargin > 3
    
    % Only plot the MAP; the ability to plot more will be added later
    for cluster = 1:nClusters
        if(iscell(clusters.assignments))
           c.selected_indices{cluster} = find(clusters.assignments{1} == cluster);
        else
           c.selected_indices{cluster} = find(clusters.assignments == cluster);
        end
    end
   
   c.projections = data.projections'; 
   c.axes = f.axes.projection_axes;
   c.acf_axis = f.axes.acf_axis;
   c.spike_times = data.spikeTimes;
   f.image = gradientClusters(f.image, c);
   p.figure{electrode} = f;
   
end