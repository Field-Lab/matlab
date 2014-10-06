function update_ica_plot
% UPDATE_ICA_PLOT   Plot the updated weights for ICADEMO
%

global signals data model prefs

%%%%%%%%%%%%%%%%%%%%%%
% plot the data points
%%%%%%%%%%%%%%%%%%%%%%
d = findobj('Tag','Original Data');

set(d,'XData',data(:,1),'YData',data(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the weight vectors
%%%%%%%%%%%%%%%%%%%%%%%%%
vect1 = findobj('Tag','Model 1 Vector');
vect2 = findobj('Tag','Model 2 Vector');

set(vect1,'XData', prefs.mean + [0 model.w1(1)], ...
	  'YData', prefs.mean + [0 model.w1(2)]);
set(vect2,'XData', prefs.mean + [0 model.w2(1)], ...
	  'YData', prefs.mean + [0 model.w2(2)]);


%%%%%%%%%%%%%%%%%%%%%%%%
% plot the input dists
%%%%%%%%%%%%%%%%%%%%%%%%

% create bin centers and format for histogram bar graph
bin_centers = [-1: 0.05: 0.95] + 0.025;

% find histogram counts and format for histogram bar graph
s1_counts = hist(signals(1,:),bin_centers);
s1_counts = s1_counts / max(s1_counts);

s2_counts = hist(signals(2,:),bin_centers);
s2_counts = s2_counts / max(s2_counts);

% grab graphic handles
s1 = findobj('Tag','Signal Distribution 1');
s2 = findobj('Tag','Signal Distribution 2');

% set the data (delete old graph and place new one)
delete(s1);
delete(s2);

% find handles for Signal axes
signal_axes = findobj('Tag','Signal Axes');
set(gcf,'CurrentAxes',signal_axes);

% create the new bar graph
hold on
s1 = bar(bin_centers, s1_counts);
s2 = bar(bin_centers, s2_counts);
set(s1,'Tag','Signal Distribution 1','FaceColor',[0 1 0]);
set(s2,'Tag','Signal Distribution 2','FaceColor',[1 0 0]);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%
% plot the output dists
%%%%%%%%%%%%%%%%%%%%%%%%

% create bin centers and format for histogram bar graph
bin_centers = [-1: 0.05: 0.95] + 0.025;

% subtract mean from data
zdata = [(data(:,1) - mean(data(:,1))) (data(:,2) - mean(data(:,2)))];

% calculate projections along each new basis vector
out1 = model.w1 * zdata';
out2 = model.w2 * zdata';

% find histogram counts and format for histogram bar graph
out1_counts = hist(out1,bin_centers);
out1_counts = out1_counts / max(out1_counts);

out2_counts = hist(out2,bin_centers);
out2_counts = out2_counts / max(out2_counts);

% grab graphic handles
out1_handle = findobj('Tag','Output Distribution 1');
out2_handle = findobj('Tag','Output Distribution 2');

% set the data (delete old graph and place new one)
delete(out1_handle);
delete(out2_handle);

% find handles for Output axes
output_axes = findobj('Tag','Output Axes');
set(gcf,'CurrentAxes',output_axes);

% create the new bar graph
hold on
out1_handle = bar(bin_centers, out1_counts);
out2_handle = bar(bin_centers, out2_counts);
set(out1_handle,'Tag','Output Distribution 1','FaceColor',[0 1 0]);
set(out2_handle,'Tag','Output Distribution 2','FaceColor',[1 0 0]);
hold off


%%%%%%%%%%%%%%%%%%%%
% update the text
%%%%%%%%%%%%%%%%%%%%

% above, calculated projections along both axes
% --> out1, out2

% calculate covariance matrix
out_covariance = cov(out1,out2);

% matrix is symmetric; either off-diagonal term is the covariance
% between out1 and out2. The diagonal terms are the auto-covariances.
out_covariance = out_covariance(2,1);

% calculate mutual information
% out_mi = mutinfo(out1,out2);

% diisplay it all
cov_text = findobj('Tag','Cov Text');
mi_text = findobj('Tag','MI Text');

% convert to text
cov_string = num2str(sprintf('%5.3g',out_covariance));
% mi_string  = num2str(sprintf('%5.3g',out_mi));

set(cov_text,'String',[' Covariance: ' cov_string]);
set(mi_text ,'String',['Mutual Info: ----- ']);
