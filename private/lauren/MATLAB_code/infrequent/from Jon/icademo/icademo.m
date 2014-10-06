
% global variables
global prefs model data

% clear the global variables
model = []; data = [];

% preferences
prefs = struct(...
    'mean',            0.5  , ...  % data: mean
    'points',         1400    ...  % data: number of points
    );

% create figure
f=figure('Position',[570 280 680 390],'MenuBar','none', ...
	 'NumberTitle','off','Name','ICA demo');

% set up GUI
x = 0.02; y = 0.83; % positioning variables


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA SELECTION CONTROL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uicontrol('Tag','Text','Style','text','String','Signal Distributions', ...
	  'Units','normalized','Position',[x y 0.2 0.05], ...
	  'BackgroundColor',[0.8 0.8 0.8]);
uicontrol('Tag','Text','Style','text','String','Signal 1', ...
          'Units','normalized','Position',[x y-0.05 0.07 0.05], ...
          'BackgroundColor',[0.8 0.8 0.8]);
uicontrol('Tag','Text','Style','text','String','Signal 2', ...
          'Units','normalized','Position',[x y-0.15 0.07 0.05], ...
          'BackgroundColor',[0.8 0.8 0.8]);

uicontrol('Tag','Signal 1 Type','Style','popupmenu','Units','normalized', ...
	  'Position', [x+0.07 y-0.05 0.13 0.06], ...
	  'String',{'Exponential','Gaussian', 'Uniform'})
uicontrol('Tag','Signal 2 Type','Style','popupmenu','Units','normalized', ...
          'Position', [x+0.07 y-0.15 0.13 0.06], ...
          'String',{'Exponential','Gaussian', 'Uniform'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AXES: SIGNAL DISTRIBUTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=y-0.40;

axes('Position',[x+0.01 y 0.2 0.2],'XTickLabel',[], ...
     'YTickLabel',[],'Color',[1 1 1],'Tag','Signal Axes');
axis([-1 1 0 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIXING MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = y - 0.11;
x = x + 0.036;
uicontrol('Tag','Text','Style','text','String','Mixing Matrix', ...
          'Units','normalized','Position',[x y+0.02 0.15 0.05], ...
          'BackgroundColor',[0.8 0.8 0.8]);
uicontrol('Tag','Text','Style','text','String','[', ...
          'Units','normalized','Position',[x-0.02 y-0.295 0.05 0.35], ...
	  'FontSize',56,'BackgroundColor',[0.8 0.8 0.8]);
uicontrol('Tag','Text','Style','text','String',']', ...
          'Units','normalized','Position',[x+0.115 y-0.295 0.05 0.35], ...
          'FontSize',56,'BackgroundColor',[0.8 0.8 0.8]);

y=y-0.05;	 
uicontrol('Tag','X_11','Style','edit','String','1', ...
	  'Units','normalized','Position',[x+0.03 y 0.035 0.05]);
uicontrol('Tag','X_12','Style','edit','String','0', ...
          'Units','normalized','Position',[x+0.085 y 0.035 0.05]);
uicontrol('Tag','X_21','Style','edit','String','0', ...
          'Units','normalized','Position',[x+0.03 y-0.1 0.035 0.05]);
uicontrol('Tag','X_22','Style','edit','String','1', ...
          'Units','normalized','Position',[x+0.085 y-0.1 0.035 0.05]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN AXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes('Position',[0.113 0.11 0.775 0.815],'XTick',[], ...
     'YTick',[],'Color',[0 0 0],'Tag','Main Axes');
axis([0 1 0 1]);
axis square;
title('ICA Demo');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL CHOICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = 0.77; y = 0.775; % positioning variables

uicontrol('Tag','Text','Style','text','String','Model Choice', ...
          'Units','normalized','Position',[x y 0.2 0.05], ...
          'BackgroundColor',[0.8 0.8 0.8]);
uicontrol('Tag','Model Choice','Style','popupmenu','Units','normalized', ...
          'Position', [x+0.02 y-0.05 0.16 0.06], 'Value', 1, ...
          'String',{'Columns of A','PCA','ICA (fixed pt)','ICA (infomax)'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AXES: OUTPUT DISTRIBUTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=y-0.35;
axes('Position',[x y 0.2 0.2],'XTickLabel',[], ...
     'YTickLabel',[],'Color',[1 1 1],'Tag','Output Axes');
axis([-1 1 0 1]);
title('Output Distributions');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT TEXT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=y-0.08;
uicontrol('Tag','Cov Text','Style','text','String',' Covariance: ------', ...
          'Units','normalized','Position',[x y 0.2 0.05], ...
          'BackgroundColor',[0.8 0.8 0.8]);
y=y-0.04;
uicontrol('Tag','MI Text','Style','text','String','Mutual Info: ------', ...
          'Units','normalized','Position',[x y 0.2 0.05], ...
          'BackgroundColor',[0.8 0.8 0.8]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE DATA BUTTON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=y-0.14;
uicontrol('Tag','Create Data','Style','pushbutton', ...
          'String','Run', ...
          'Units','normalized','Position',[x+0.05 y 0.1 0.08], ...
          'Callback','generate_ica_data');





% Start the plots
init_ica_plot;
