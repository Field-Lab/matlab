function gui_cb (action)

% This file is used by FASTICAG

% This file holds the callbacks to the main window

% 3.4.1998
% Hugo Gävert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables

% Handle to the main figure
global h_f_FastICA;

% Handles for needed controls in main figure;
global h_t_mixedStatus;
global h_t_dim;
global h_t_numOfSamp;
global h_t_newDim;
global h_t_whiteStatus;
global h_t_icaStatus;
global h_pm_approach;
global h_e_numOfIC;
global h_pm_g;

% Global variable used by gui_cb & gui_lc
global g_mix_meanremoved;

% Main values are stored here
global g_FastICA_mixedsig;
global g_FastICA_pca_D;
global g_FastICA_pca_E;
global g_FastICA_white_sig;
global g_FastICA_white_wm;
global g_FastICA_white_dwm;
global g_FastICA_ica_sig;
global g_FastICA_ica_A;
global g_FastICA_ica_W;
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

% String values are here
global c_FastICA_appr_strV;
global c_FastICA_g1_strD;
global c_FastICA_g1_strV;
global c_FastICA_g2_strD;
global c_FastICA_g2_strV;
global c_FastICA_iSta_strV;
global c_FastICA_dMod_strV;
global c_FastICA_verb_strV;

% What is the load type of load dialog
global g_FastICA_loadType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What ever we do, it will take some time... not much, but some :-)
watchonInFigure = watchon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(action, 'InitAll')

  % If the data is already loaded, then get the information from data
  % and show to the user (also set g_FastICA_numOfIC)
  if ~isempty(g_FastICA_mixedsig)
    set(h_t_mixedStatus, 'String', '');
    [dim, numofsamp] = size(g_FastICA_mixedsig);
    set(h_t_dim, 'String', int2str(dim));
    set(h_t_numOfSamp, 'String', int2str(numofsamp));
    set(h_t_newDim, 'String', int2str(dim));
    set(h_e_numOfIC, 'String', int2str(dim));
    g_FastICA_numOfIC = dim;
    g_mix_meanremoved = 0;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'LoadData')

  handle = findobj('Tag','f_FastICALoad');  % Check if the window is already
  if isempty(handle)                        % open. If not then open it.
    pos = get(h_f_FastICA, 'Position');
    gui_l(pos(1), pos(2));
  else
    if strcmp(g_FastICA_loadType, 'data')   % Check if it was the same load
      figure(handle);                       % window. If it wasn't then
    else                                    % close the other window first
      close(handle);                        % and then open the load window
      fprintf('''Load initial guess'' -dialog closed!\n');
      pos = get(h_f_FastICA, 'Position');
      gui_l(pos(1), pos(2));
    end
  end

  % gui_cb NewData; - is called from the load function if not canceled...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'NewData')

  % New data is loaded or the old data changed. We need to find out
  % somethings about the new data... and do some other stuff also...
  [dim, numofsamp] = size(g_FastICA_mixedsig);
  set(h_t_dim, 'String', dim);
  set(h_t_newDim, 'String', dim);
  set(h_t_numOfSamp, 'String', numofsamp);
  set(h_e_numOfIC, 'String', int2str(dim));

  g_FastICA_numOfIC = dim;    % Default for numOfIC = the new dimension
                              % PCA needs to be calculated again.
  g_FastICA_pca_E = [];       % We use this to check if PCA is calculated
  gui_cb NullWhite;           % Whitening needs to be done again also

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'NullWhite')

                              % Whitening needs to done again next time it's needed
  g_FastICA_white_sig = [];   % We use this to check if whitening is calculated
  gui_cb NullICA;             % The FPICA must be calculated again

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'NullICA')

                              % If IC's are needed they have to bee calculated again
  g_FastICA_ica_sig = [];     % We use this to check if FPICA is calculated
  set(h_t_icaStatus,'String','Not yet done');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'Transpose')

  g_FastICA_mixedsig = g_FastICA_mixedsig';
  gui_cb NewData;             % Data has been changed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'DoPCA')

  if ~isempty(g_FastICA_mixedsig)
    % We'll remove mean of the data here also just in case...
    if ~(g_mix_meanremoved)
      g_mix_meanremoved = 1;
      g_FastICA_mixedsig = remmean(g_FastICA_mixedsig);
    end

    % Do PCA interactively: ask the user for eigenvalues
    [g_FastICA_pca_E, g_FastICA_pca_D] = pcamat(g_FastICA_mixedsig, ...
          0, 0, 'gui', deblank(c_FastICA_verb_strV(g_FastICA_verbose,:)));

    newdim = size(g_FastICA_pca_D, 1);
    set(h_t_newDim, 'String', int2str(newdim));
    set(h_e_numOfIC, 'String', int2str(newdim));
    g_FastICA_numOfIC = newdim;
    gui_cb NullWhite;           % Whitening needs to be done again also
                                % but we'll do it when it's needed.
  else
    fprintf('Data not loaded yet!\n\n');
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'OrigDim')

  gui_cb NewData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'ShowMixed')

  if ~isempty(g_FastICA_mixedsig)
    handle = findobj('Tag','f_FastICA_mix');  % Check if the window is already
    if isempty(handle)                        % open. If not then open it.
      figure('Tag', 'f_FastICA_mix', ...
             'Name', 'FastICA: Plot data', ...
             'NumberTitle', 'off');
    else
      figure(handle);
      clf;		% clear the figure for next plots
    end
    dispsig(g_FastICA_mixedsig');
    if g_mix_meanremoved
      title('Mixed signals (mean removed)');
    else
      title('Mixed signals');
    end
  else
    fprintf('Data not loaded yet!\n\n');
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'ShowWhite')

  if ~isempty(g_FastICA_mixedsig)
    if isempty(g_FastICA_white_sig)     % if whitening is not done, we need to
      gui_cb Whiten;                    % do it before we can display the
    end                                 % whitened signals

    handle = findobj('Tag','f_FastICA_white');  % Check if the window is already
    if isempty(handle)                          % open. If not then open it.
      figure('Tag', 'f_FastICA_white', ...
             'Name', 'FastICA: Plot whitened', ...
             'NumberTitle', 'off');
    else
      figure(handle);
      clf;		% clear the figure for next plots
    end
    dispsig(g_FastICA_white_sig');
    title('Whitened signals');
  else
    fprintf('Data not loaded yet!\n\n');
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp (action, 'Whiten')

  set(h_t_whiteStatus,'String','Computing...');

  % If PCA is not calculated, we'll have to calculate it now,
  % we'll do it without guestions - we don't reduce the dimension
  % here - but PCAMAT might reduce the dimension automatically.
  if isempty(g_FastICA_pca_E)
    % We'll remove mean of the data here also just in case...
    if ~(g_mix_meanremoved)
      g_mix_meanremoved = 1;
      g_FastICA_mixedsig = remmean(g_FastICA_mixedsig);
    end

    [g_FastICA_pca_E, g_FastICA_pca_D] = pcamat(g_FastICA_mixedsig, ...
                     1, size(g_FastICA_mixedsig, 1), 'off', ...
                     deblank(c_FastICA_verb_strV(g_FastICA_verbose,:)));

    % Check if the dimension was reduced automatically
    newdim = size(g_FastICA_pca_D, 1);
    set(h_t_newDim, 'String', int2str(newdim));
    % Check if the numOfIC now has illegal value entered 
    % We do that by telling the program that there is new value 
    % entered for NumOfIC.
    gui_cb ChangeNumOfIC;
  end

  % And now we can calculate whitening...
  [g_FastICA_white_sig, g_FastICA_white_wm, g_FastICA_white_dwm] = ...
     whitenv(g_FastICA_mixedsig, g_FastICA_pca_E, g_FastICA_pca_D, ...
             deblank(c_FastICA_verb_strV(g_FastICA_verbose,:)));

  set (h_t_whiteStatus,'String','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'ChangeApproach')

  % Get the old value for g
  eval(['g_str = c_FastICA_g' int2str(g_FastICA_approach) '_strV;']);
  old_g = deblank(g_str(g_FastICA_g,:));

  % Get and set the new value for approach
  g_FastICA_approach = get(h_pm_approach, 'Value');

  % The possible values for g depend on the value of approach...
  eval(['g_str = c_FastICA_g' int2str(g_FastICA_approach) '_strD;']);
  set(h_pm_g, 'String', g_str);

  % Match the old g value from the new g values so that if the 
  % old_g can be found from the new values (anywhere), then set new g
  % to that value, and if it's not found then set the new value to 1.
  match = 0;
  eval(['g_str = c_FastICA_g' int2str(g_FastICA_approach) '_strV;']);
  for i=1:size(g_str,1)
    if strcmp(old_g, deblank(g_str(i,:)))
      match = i;
    end
  end
  if match == 0
    match = 1;   % the old g is not availabe anymore, set g = 1.
  end
  g_FastICA_g = match;
  set(h_pm_g, 'Value', match);

  % We must also check if the numOfIC can be awailable.
  % If symmetric...
  if strcmp(deblank(c_FastICA_appr_strV(g_FastICA_approach,:)),'symm')
    % The numOfIC must be equal to newDim
    set(h_e_numOfIC, 'String', get(h_t_newDim, 'String'), 'Enable', 'off');
    g_FastICA_numOfIC = str2num(get(h_t_newDim, 'String'));

  % else deflation...
  else
    % The numOfIC can be different than newDim
    set(h_e_numOfIC, 'Enable', 'on');
  end

  gui_cb NullICA;      % The options are changed so we must calculate ICA again

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'ChangeNumOfIC')

  % Get the new value... and store it later on after some checks
  numofic = str2num(get(h_e_numOfIC, 'String'));

  % The number of IC can't be less than 1 or more than the reduced dimension.
  numoficmax = str2num(get(h_t_newDim, 'String'));
  if numofic < 1
    set(h_e_numOfIC, 'String', '1');
    g_FastICA_numOfIC = 1;
  elseif numofic > numoficmax
    set(h_e_numOfIC, 'String', int2str (numoficmax));
    g_FastICA_numOfIC = numoficmax;
  else
    g_FastICA_numOfIC = numofic;
  end

  gui_cb NullICA;      % The options are changed so we must calculate ICA again

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'ChangeG')

  % Get the new value for g.
  g_FastICA_g = get(h_pm_g, 'Value');

  gui_cb NullICA;      % The options are changed so we must calculate ICA again

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'AdvOpt')

  handle = findobj('Tag','f_FastICAAdvOpt');
  if isempty(handle)                        % Check to see if the window is
    pos = get(h_f_FastICA, 'Position');     % already open...
    gui_adv(pos(1), pos(2));
  else
    figure(handle)
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'ShowICASig')

  if ~isempty(g_FastICA_mixedsig)
    handle = findobj('Tag','f_FastICA_ica');  % Check if the window is already
    if isempty(handle)                        % open. If not then open it.
      figure('Tag', 'f_FastICA_ica', ...
             'Name', 'FastICA: Plot ICs', ...
             'NumberTitle', 'off');
    else
      figure(handle);
      clf;		% clear the figure for next plots
    end

    % If the IC's are not already calculated, we'll do it now
    if isempty(g_FastICA_ica_sig)
      gui_cb DoFPICA;

      % If the value for DisplayMode was not 'none' then the signals
      % are allready displayed. 
      if strcmp(deblank(c_FastICA_dMod_strV(g_FastICA_displayMo,:)),'off')
        dispsig(g_FastICA_ica_sig')
      end
    else
      dispsig(g_FastICA_ica_sig');      % The FPICA was already done, so we just
    end                                 % display the signals
    title('Independent components');
  else
    fprintf('Data not loaded yet!\n\n');
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(action, 'DoFPICA')

  if ~isempty(g_FastICA_mixedsig)
    if isempty(g_FastICA_white_sig)     % We need the whitened signal here
      gui_cb Whiten;
    end

    set(h_t_icaStatus,'String','Computing...');

    % The possible values for g depend on approach
    eval(['g_str = c_FastICA_g' int2str(g_FastICA_approach) '_strV;']);

    % We'll contruct a command string which we'll later on evaluate...
    % This is where the Fixed point algorithm is used.
    command_str = ['[g_FastICA_ica_A,g_FastICA_ica_W]=' ...
      'fpica(g_FastICA_white_sig, g_FastICA_white_wm, g_FastICA_white_dwm,' ...
      '''' deblank(c_FastICA_appr_strV(g_FastICA_approach,:)) ''',' ...
      'g_FastICA_numOfIC,' ...
      '''' deblank(g_str(g_FastICA_g,:)) ''',' ...
      'g_FastICA_a1, g_FastICA_a2, g_FastICA_epsilon, g_FastICA_maxNumIte,' ...
      '''' deblank(c_FastICA_iSta_strV(g_FastICA_initState,:)) ''',' ...
      'g_FastICA_initGuess,' ...
      '''' deblank(c_FastICA_dMod_strV(g_FastICA_displayMo,:)) ''',' ...
      'g_FastICA_displayIn,' ...
      '''' deblank(c_FastICA_verb_strV(g_FastICA_verbose,:)) ''');'];

    % If the user wants to plot while computing...
    % let's at least plot it to the right figure then
    if ~strcmp(deblank(c_FastICA_dMod_strV(g_FastICA_displayMo,:)),'off')
      handle = findobj('Tag','f_FastICA_ica');  % Check if the window is already
      if isempty(handle)                        % open. If not then open it.
        figure('Tag', 'f_FastICA_ica', ...
               'Name', 'FastICA: Plot ICs', ...
               'NumberTitle', 'off');
      else
        figure(handle);
        clf;		% clear the figure for next plots
      end
    end

    % ... and so let's do it...
    eval(command_str);

    g_FastICA_ica_sig = g_FastICA_ica_W * g_FastICA_mixedsig;

    set (h_t_icaStatus,'String','Done');
  else
    fprintf('Data not loaded yet!\n\n');
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp (action, 'SaveData')

  handle = findobj('Tag','f_FastICASave');  % Check if the window is already
  if isempty(handle)                        % open. If not then open it.
    pos = get(h_f_FastICA, 'Position');
    gui_s(pos(1), pos(2));
  else
    figure(handle);                       % window. If it wasn't then
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp (action, 'Quit')

  % We'll close the other dialogs if they are open.
  Tags = ['f_FastICALoad  '
          'f_FastICAAdvOpt'
          'f_FastICASave  '
          'f_FastICA_mix  '
          'f_FastICA_white'
          'f_FastICA_ica  '];
  for i=1:size(Tags,1)
    handle = findobj('Tag', deblank(Tags(i,:)));
    if ~isempty(handle)
      close(handle);
    end
  end

  % Close this window
  close(h_f_FastICA);

  % Clear the used global variables.
  gui_cg;

  % Use break to 'jump' over the watchoff statement at the end
  break;
  % ... and we're done.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp (action, 'About')

  gui_help('gui_cb_about');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp (action, 'Help')

  gui_help('gui_cb_help');
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end    % if ... elseif ...

watchoff (watchonInFigure);
