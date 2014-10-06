function gui_lc (action)

% This file is used by FASTICAG

% This file holds the callbacks for load-dialog

% 3.4.1998
% Hugo Gävert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables

% Handle to the window
global h_f_FastICALoad;

% Global variable used by gui_cb & gui_lc
global g_mix_meanremoved;

% Handles to some of the controls in window
global h_e_file;

% Needed handles from the main figure
global h_t_mixedStatus;
% Needed handles from the advOpt figure
global h_b_initGuess;
global h_t_initGuess;
global h_pm_initState;

% The needed main variables
global g_FastICA_mixedsig;

% What is the load type of load dialog
global g_FastICA_loadType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This should not take long...
watchonInFigure = watchon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(action, 'Load')

  varName = get(h_e_file, 'String');      % The name of the variable to be loaded
  command=['evalin(''base'',''assignin(''''caller'''',''''data'''',' varName ')'')'];
  eval(command,'fprintf(''Variable not found in MATLAB workspace, data not loaded!\n'');data=[];');                          % Variable is copyed to 'data'

  if isempty(data)  % if there was no name given...
    watchoff (watchonInFigure);
    break;
  end
  if strcmp(g_FastICA_loadType,'data')          % New data
    g_FastICA_mixedsig = data;
    set(h_t_mixedStatus, 'String', '');
    g_mix_meanremoved=0;                        % New data - so that means ...
    gui_cb NewData;                             

  elseif strcmp(g_FastICA_loadType,'guess')     % New initial guess
    set(h_b_initGuess, 'UserData', data);       % Since we loaded new initial
    set(h_t_initGuess, 'String', 'Loaded');     % guess, we wan't to use it too
    set(h_pm_initState, 'Value', 2);		% ... set initState to 'guess'

  end

  close(h_f_FastICALoad);                       % close the dialog

  % Use break to 'jump' over the watchoff statement at the end
  break;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp (action, 'Cancel')
  close(h_f_FastICALoad);                       % do nothing just exit

  % Use break to 'jump' over the watchoff statement at the end
  break;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'Help')

  % Which help do we show?
  if strcmp(g_FastICA_loadType,'data')
    gui_help('gui_lc_data');
  elseif strcmp(g_FastICA_loadType,'guess')   
    gui_help('gui_lc_guess');
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end    % if ... elseif ...

watchoff (watchonInFigure);
