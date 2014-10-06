function init_ica_plot
% INIT_ICA_PLOT    Initializes graphic handles and axes for ICADEMO
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN AXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab the handle for the plot
axes_handle = findobj('Tag','Main Axes');
set(gcf,'CurrentAxes',axes_handle);

hold on;

% -- Data points
plot(-1,-1,'b.','MarkerSize',6,'Tag','Original Data');

% -- Weight vectors
plot(-1,-1,'g-','LineWidth',4,'Tag','Model 1 Vector');
plot(-1,-1,'r-','LineWidth',4,'Tag','Model 2 Vector');

hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNAL DISTRIBUTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab the handle for the plot
axes_handle = findobj('Tag','Signal Axes');
set(gcf,'CurrentAxes',axes_handle);

hold on;

% -- plot histograms
s1 = bar(0,0);
s2 = bar(0,0);
set(s1, 'Tag','Signal Distribution 1');
set(s2, 'Tag','Signal Distribution 2');

hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT DISTRIBUTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab the handle for the plot
axes_handle = findobj('Tag','Output Axes');
set(gcf,'CurrentAxes',axes_handle);


hold on;

% -- plot histograms
o1 = bar(0,0);
o2 = bar(0,0);
set(o1, 'Tag','Output Distribution 1');
set(o2, 'Tag','Ourput Distribution 2');

hold off;



