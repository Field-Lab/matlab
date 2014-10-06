function generate_ica_data
% GENERATE_ICA_DATA:  Generates data for ICADEMO
%

global data prefs model signals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. CLEAR GLOBAL VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data =[]; model = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. CREATE THE SIGNALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab GUI choices 
s1_menu = findobj('Tag','Signal 1 Type');
s2_menu = findobj('Tag','Signal 2 Type');

s1_type = get(s1_menu,'Value');
s2_type = get(s2_menu,'Value');

% create the signals
signals = zeros(2,prefs.points);

signals(1,:) = generate_ica_dist(s1_type,prefs.points);
signals(2,:) = generate_ica_dist(s2_type,prefs.points);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. MIX THE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grab the mixing maxtrix
A = mix_matrix;

% Multiply by mixing matrix
xy = A * signals;

% add the mean
xy = xy' + prefs.mean;

% set the data
data = xy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. UPDATE MODEL & PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% update the model (and plot)
generate_ica_model;


