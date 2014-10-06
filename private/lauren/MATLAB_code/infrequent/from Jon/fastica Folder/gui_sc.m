function gui_sc (action)

% This file is used by FASTICAG

% This file holds the callbacks for save-dialog

% 1.4.1998
% Hugo Gävert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables

% Handle to the window
global h_f_FastICASave;

% Handles to some of the controls in window
global h_e_suffix;

% The needed main variables
global g_FastICA_ica_sig;
global g_FastICA_ica_A;
global g_FastICA_ica_W;
global g_FastICA_white_sig;
global g_FastICA_white_wm;
global g_FastICA_white_dwm;
global g_FastICA_pca_E;
global g_FastICA_pca_D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This should not take long...
watchonInFigure = watchon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(action, 'Save')

  suffix = deblank(get(h_e_suffix, 'String')); % The suffix for the variables
%  if ~isempty(suffix)
%    suffix = ['_' suffix];
%  end

  fprintf('Saving results in variables in Matlab workspace.\n');
  assignin('base',['IC' suffix],g_FastICA_ica_sig);
  assignin('base',['A' suffix],g_FastICA_ica_A);
  assignin('base',['W' suffix],g_FastICA_ica_W);
  assignin('base',['whitesig' suffix],g_FastICA_white_sig);
  assignin('base',['whiteningMatrix' suffix],g_FastICA_white_wm);
  assignin('base',['dewhiteningMatrix' suffix],g_FastICA_white_dwm);
  assignin('base',['E' suffix],g_FastICA_pca_E);
  assignin('base',['D' suffix],g_FastICA_pca_D);

  close(h_f_FastICASave);                  % close the dialog

  % Use break to 'jump' over the watchoff statement at the end
  break;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp (action, 'Cancel')
  close(h_f_FastICASave);                       % do nothing just exit

  % Use break to 'jump' over the watchoff statement at the end
  break;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'Help')

  gui_help('gui_sc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end    % if ... elseif ...

watchoff (watchonInFigure);
