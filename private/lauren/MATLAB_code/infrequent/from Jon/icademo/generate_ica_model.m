function generate_ica_model
% GENERATE_ICA_MODEL:  Initializes model in ICADEMO
%

global data model

% set up model data structure
model = struct('w1',[],'w2',[]);

% grab selected learning rule
popupmenu = findobj('Tag','Model Choice');
menu_item = get(popupmenu,'Value');

switch menu_item
 case 1
  % Columns of A
  % ============

  % grab the mixing matrix
  A = mix_matrix;

  % find normalization (to make it pretty)
  a_max = max(max(A)) / 0.3;

  % set the model to the mixing matrix
  model.w1 = A(:,1)' / a_max;
  model.w2 = A(:,2)' / a_max;


 case 2
  % PCA
  % ===============

  % remove mean
  mn = mean(data);
  newdata(:,1) = data(:,1) - mn(1);
  newdata(:,2) = data(:,2) - mn(2);

  % compute svd
  [u,s,v] = svd(newdata);
  % "normalize" according to variance
  s = diag(s); s = s / s(1) / 2;

  % set weights
  model.w1 = v(1,:) * s(1);
  model.w2 = v(2,:) * s(2);

 case 3
  % ICA (fixed point)
  % =================

  [junk,A,W] = fastica(data','g','tanh','approach','symm','epsilon', ...
		       0.0001,'displayMode','off','verbose','off');

  % set weights
  model.w1 = A(:,1)';
  model.w2 = A(:,2)';

 case 4
  % ICA (infomax)
  % ==============
  
  [junk,A,W] = infomax(data');

  % set weights
  model.w1 = A(:,1)';
  model.w2 = A(:,2)';

end

% update the plot
update_ica_plot;

