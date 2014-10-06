function gui_advc (action)

% This file is needed by FASTICAG

% This file holds the callbacks for advanced options -dialog

% 1.4.1998
% Hugo Gävert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables

% Handle to the window
global h_f_FastICAAdvOpt;

% Handles to some of the controls in window
global h_e_a1;
global h_e_a2;
global h_e_epsilon;
global h_e_maxIterations;
global h_pm_initState;
global h_b_initGuess;
global h_t_initGuess;
global h_pm_displayMode;
global h_e_displayInterval;
global h_pm_verbose;

% Needed handles to the main window
global h_f_FastICA;
global h_t_icaStatus;

% Some of the main variables needed
global g_FastICA_ica_sig;
global g_FastICA_initGuess;
global g_FastICA_approach;
global g_FastICA_numOfIC;
global g_FastICA_g;
global g_FastICA_a1;
global g_FastICA_a2;
global g_FastICA_epsilon;
global g_FastICA_maxNumIte;
global g_FastICA_initState;
global g_FastICA_displayMo;
global g_FastICA_displayIn;
global g_FastICA_verbose;

global c_FastICA_iSta_strV;

% What is the load type of load dialog
global g_FastICA_loadType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This should not take long...
watchonInFigure = watchon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(action, 'Checka1')
  e_a1_val = str2num(get(h_e_a1, 'String'));
  if e_a1_val <= 0
    set(h_e_a1, 'String', '0.1');
%  elseif e_a1_val > 1
%    set(h_e_a1, 'String', '1');
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'Checka2')
  e_a2_val = str2num(get(h_e_a2, 'String'));
  if e_a2_val <= 0
    set(h_e_a2, 'String', '0.1');
%  elseif e_a2_val > 1
%    set(h_e_a2, 'String', '1');
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'loadGuess')

  handle = findobj('Tag','f_FastICALoad');  % Check if the window is already
  if isempty(handle)                        % open. If not then open it.
    pos = get(h_f_FastICA, 'Position');
    gui_l(pos(1), pos(2));
  else
    if strcmp(g_FastICA_loadType, 'guess')  % Check if it was the same load
      figure(handle);                       % window. If it wasn't then
    else                                    % close the other window first
      close(handle);                        % and then open the load window
      fprintf('''Load data'' -dialog closed!\n');
      pos = get(h_f_FastICA, 'Position');
      gui_l(pos(1), pos(2));
    end
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'OK')
  newValues = 0;

  val = g_FastICA_a1;
  g_FastICA_a1 = str2num(get(h_e_a1, 'String'));
  if (g_FastICA_a1 ~= val)
    newValues = 1;
  end

  val = g_FastICA_a2;
  g_FastICA_a2 = str2num(get(h_e_a2, 'String'));
  if (g_FastICA_a2 ~= val)
    newValues = 1;
  end

  val = g_FastICA_epsilon;
  g_FastICA_epsilon = str2num(get(h_e_epsilon, 'String'));
  if (g_FastICA_epsilon ~= val)
    newValues = 1;
  end

  val = g_FastICA_maxNumIte;
  g_FastICA_maxNumIte = str2num(get(h_e_maxIterations, 'String'));
  if (g_FastICA_maxNumIte ~= val)
    newValues = 1;
  end

  val = g_FastICA_initState;
  g_FastICA_initState = get(h_pm_initState, 'Value');
  if (g_FastICA_initState ~= val)
    newValues = 1;
  end

  val = g_FastICA_initGuess;
  g_FastICA_initGuess = get(h_b_initGuess, 'UserData');
  if min(size(val) == size(g_FastICA_initGuess)) == 0
    newValues = 1;
  else
    if (g_FastICA_initGuess ~= val)
      newValues = 1;
    end
  end

  g_FastICA_displayMo = get(h_pm_displayMode, 'Value');
  g_FastICA_displayIn = str2num(get(h_e_displayInterval, 'String'));
  g_FastICA_verbose = get(h_pm_verbose, 'Value');

  if newValues == 1
    gui_cb NullICA;
  end
  close(h_f_FastICAAdvOpt);

  % Use break to 'jump' over the watchoff statement at the end
  break;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp (action, 'Cancel')
  close(h_f_FastICAAdvOpt);

  % Use break to 'jump' over the watchoff statement at the end
  break;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'Default')
  % set default values to controls
  set(h_e_a1, 'String','1');
  set(h_e_a2, 'String','1');
  set(h_e_epsilon, 'String','0.0001');
  set(h_e_maxIterations, 'String','1000');
  set(h_pm_initState, 'Value',1);
  set(h_pm_displayMode, 'Value',1);
  set(h_e_displayInterval, 'String','1');
  set(h_pm_verbose, 'Value',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'Help')

  gui_help('gui_advc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end    % if ... elseif ...

watchoff (watchonInFigure);