function gui_help(which_help)

% Used by FASTICAG

% All the help texts and title used by GUI are stored here. 
% Make changes here.
% Also displays the helpwindow with the selected text

% Hugo Gävert
% 23.4.1998



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(which_help, 'pcamat')

  helptitle = 'FastICA: Reduce dimension';
  helptext=[ ...
    'You may reduce the dimension of the data by selecting only the     '
    'subspace corresponding to certain eigenvalues of the covariance    '
    'matrix of the data. Give the indices of the first and last         '
    'eigenvalues (sorted in descending order) to be included (all       '
    'eigenvalues in between will be included as well).  The eigenvalues '
    'and their indices can be seen in the graphical plot now on the     '
    'screen. The heights of the bars give the eigenvalues, with indices '
    'below.                                                             '
    '                                                                   '
    'For example, give ''1'' and ''n'' if you want to reduce the dimension  '
    'to n by principal component analysis, which means discarding the   '
    'subspaces corresponding to the smallest eigenvalues. Such a        '
    'dimension reduction may reduce noise and improve the performance of'
    'ICA.                                                               '];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(which_help, 'gui_cb_about')

  helptitle='About FastICA';
  helptext =[ ...
    'FastICA for Matlab 5.x                                                    '
    'Version 1.01  April 23 1998                                               '
    'Copyright (c) Jarmo Hurri, Hugo Gävert, Jaakko Särelä, and Aapo Hyvärinen.'
    '                                                                          '
    'For more information please see:                                          '
    'http://www.cis.hut.fi/projects/ica/fastica/                               '];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(which_help, 'gui_cb_help')

  helptitle='FastICA GUI';
  helptext = [...
    'Basic function:                                                            '
    '                                                                           '
    '- Click LOAD DATA and give the name of the variable that contains          '
    '  the data.                                                                '
    '                                                                           '
    '- Click DO ICA to perform the analysis.                                    '
    '                                                                           '
    '- Click SAVE RESULTS to store the results for future use.                  '
    '                                                                           '
    'Options:                                                                   '
    '                                                                           '
    'If the input matrix contains the signals as column vectors instead of      '
    'row vectors, click on TRANSPOSE to transpose the data matrix.              '
    '                                                                           '
    'Click on PLOT DATA to see the data as 1-D time signals.                    '
    '                                                                           '
    'Clicking REDUCE DIM gives you a graphical plot of the eigenvalue           '
    'structure of the covariance matrix of the data. You can then reduce        '
    'the dimension of the data by retaining only the subspaces corresponding to '
    'the largest (or smallest) eigenvalues (i.e. variances). To undo this       '
    'operation click ORIGINAL DIM. You can plot the whitened (preprocessed      '
    'data) by PLOT WHITENED.                                                    '
    '                                                                           '
    'Click on DO ICA to perform independent component analysis.                 '
    'Clicking on PLOT ICS has the same effect, except that DO ICA forces        '
    'recomputation of ICA.                                                      '
    '                                                                           '
    'You can choose the decorrelation approach by the ''approach'' drop-down menu:'
    'deflation means that the independent components are estimated              '
    'one-by-one, whereas in the symmetric approach they are estimated in        '
    'parallel. In the deflation approach, you can also choose the number of     '
    'independent component to be estimated (it must not be larger than the      '
    'dimension of the data).                                                    '
    '                                                                           '
    'You have a choice of three nonlinearities:                                 '
    '                                                                           '
    '''pow3'' (default) :  g(u)=u^3                                               '
    '''tanh''           :  g(u)=tanh(u)                                           '
    '''gauss''          :  g(u)=u*exp(-u^2/2)                                     '
    '                                                                           '
    'For example, you could choose approach=''symmetric'' and nonlinearity=''tanh'' '
    'to perform maximum likelihood ICA estimation for supergaussian data.       '
    '                                                                           '
    'The ADVANCED OPTIONS menu has its own HELP button.                         '
    '                                                                           '
    'To interrupt FastICA during computations, press CTRL+C in the main         '
    'MATLAB window. Use this with caution, since it stops the execution of      '
    'the main program as well. Therefore the main program may get confused      '
    '(for example, the mouse pointer will appear incorrectly).                  '];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(which_help, 'gui_advc')

  helptitle='FastICA GUI: Advanced options';
  helptext = [...
    'Advanced options:                                                     '
    '                                                                      '
    'You can fine-tune the nonlinearity used in the fixed-point algorithm. '
    'The nonlinearities tanh and gauss contain parameters a1 and a2, so    '
    'that the nonlinearities are in fact defined as:                       '
    '''tanh''           :  g(u)=tanh(a1*u)                                   '
    '''gauss''          :  g(u)=u*exp(-a2*u^2/2)                             '
    'The default values of a1 and a2 are 1, in which case they effectively '
    'disappear from the definitions.                                       '
    '                                                                      '
    'The parameter ''epsilon'' is used to decide if the algorithm has        '
    'converged. A larger epsilon makes the convergence test less strict.   '
    '                                                                      '
    '''Maximum number of iterations'' gives the absolute maximum of          '
    'iterations used in the estimation procedure. In the deflation         '
    'approach, this is iterations per component.                           '
    '                                                                      '
    'You can input the ''Initial state'' of the algorithm, i.e. the initial  '
    'value for A. Choose ''guess'' in the drop-down menu ''Initial state'',    '
    'click on ''Load Initial guess'', and give the name of the variable in   '
    'Matlab workspace that contains the initial value. Only possible when  '
    '''approach'' is ''symmetric'' in the main menu.                           '
    '                                                                      '
    'In the drop-down menu ''Display mode'' you can choose if the results are'
    'plotted.  You may wish to switch this off especially if you have lots '
    'of data which takes a long time to plot.  ''Iteration between displays'''
    'tells how often the running estimates of the independent components   '
    'are plotted: A value of 1 means after every iteration.                '
    '                                                                      '
    'Click on DEFAULT to return to default values for all advanced options.'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(which_help, 'gui_lc_data')
  helptitle='FastICA GUI: Load data';
  helptext = [...
    'Input the name of the variable in Matlab workspace that contains the  '
    'data. The data must be in a single matrix, each row (or column) giving'
    'the values of one signal. If the signals are in column vectors, click '
    'TRANSPOSE after loading the data to transpose the data matrix.        '
    '                                                                      '
    'If the data is in a file, load it to Matlab workspace first.          '];

elseif strcmp(which_help, 'gui_lc_guess')
  helptitle='FastICA GUI: Load guess';
  helptext = [...
    'Input the name of the variable in Matlab workspace that contains the'
    'initial value for the mixing matrix A, and click OK. If the initial '
    'value is in a file, load it to Matlab workspace first.              '];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(which_help, 'gui_sc')

  helptitle='FastICA GUI: Save results';
  helptext = [...
    'The results will be saved as variables in Matlab workspace.            '
    'You give a suffix that identifies the variables. For example, if you   '
    'give ''_FASTICA'', the results will be stored in the following variables:'
    '                                                                       '
    'W_FASTICA   : estimate of the separating matrix                        '
    'A_FASTICA   : estimate of the mixing matrix, obtained as the pseudo-   '
    '              inverse of the separating matrix                         '
    'IC_FASTICA  : estimated independent components (row vectors)           '
    '                                                                       '
    'Additional results related to preprocessing:                           '
    'D_FASTICA and E_FASTICA    : give the eigenvalue decomposition of the  '
    '                             covariance matrix                         '
    'whiteningMatrix_FASTICA    : matrix performing whitening and dimension '
    '                             reduction                                 '
    'dewhiteningMatrix_FASTICA  : the pseudoinverse of the whitening matrix '
    'whitesig_FASTICA           : whitened (i.e. preprocessed) signals.     '];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

helpwin(helptext, helptitle);