function [Out1, Out2, Out3] = fastica(mixedsig, s1, v1, s2, v2, s3, v3, ...
      s4, v4, s5, v5, s6, v6, s7, v7, s8, v8, ...
      s9, v9, s10, v10, s11, v11, s12, v12, ...
      s13, v13, s14, v14, s15, v15, s16, v16, ...
      s17, v17, s18, v18, s19, v19, s20, v20);

% FASTICA - Fast Independent Component Analysis.
%
% FASTICA(mixedsig) estimates the independent components from given
% multidimensional signals. Each row of matrix mixedsig is one
% observed signal.  FASTICA uses Hyvarinen's fixed-point algorithm,
% see http://www.cis.hut.fi/projects/ica/fastica/. Output from the
% function depends on the number output arguments:
%
% [icasig] = FASTICA (mixedsig); the rows of icasig contain the
% estimated independent components.
%
% [icasig, A, W] = FASTICA (mixedsig); outputs the estimated separating
% matrix W and the corresponding mixing matrix A (obtained as
% pseudoinverse of W) as well.
%
% [A, W] = FASTICA (mixedsig); gives only the estimated mixing matrix
% A and the separating matrix W.
%
% Some optional arguments induce other output formats, see below.
%
% A graphical user interface for FASTICA can be launched by the
% command FASTICAG
%
% FASTICA can be called with numerous optional arguments. Optional
% arguments are given in parameter pairs, so that first argument is
% the name of the parameter and the next argument is the value for
% that parameter. Optional parameter pairs can be given in any order.
%
% OPTIONAL PARAMETERS:
%
% Parameter name        Values and description
%
%======================================================================
% --Basic parameters in fixed-point algorithm:
%
% 'approach'            (string) The decorrelation approach used. Can be
%                       symmetric ('symm'), i.e. estimate all the
%                       independent component in parallel, or
%                       deflation ('defl'), i.e. estimate independent
%                       component one-by-one like in projection pursuit.
%                       Default is 'defl'.
%
% 'numOfIC'             (integer) Number of independent components to
%                       be estimated.  Ignored if 'approach'='symm'.
%                       Default equals the dimension of data.
%
%======================================================================
% --Choosing the nonlinearity:
%
% 'g'                   (string) Chooses nonlinearity g used in 
%                       the fixed-point algorithm. Possible values:
%
%                       Value of 'g':      Nonlinearity used:
%                       'pow3' (default)   g(u)=u^3
%                       'tanh'             g(u)=tanh(a1*u)
%                       'gaus'             g(u)=u*exp(-a2*u^2/2)
%
% 'a1'                  (number) Parameter a1 used when g='tanh'.
%                       Default is 1.
% 'a2'                  (number) Parameter a2 used when g='gaus'.
%                       Default is 1.
% 
%======================================================================
% --Controlling convergence:
%
% 'epsilon'             (number) Stopping criterion. Default is 0.0001.
%
% 'maxNumIterations'    (integer) Maximum number of iterations.
%                       Default is 1000.
%
% 'initGuess'           (matrix) Initial guess for A. Default is random.
%
%======================================================================
% --Graphics and text output:
%
% 'verbose'             (string) Either 'on' or 'off'. Default is
%                       'on': report progress of algorithm in text format.
%                     
% 'displayMode'         (string) Plot running estimates of independent
%                       components: 'on' or 'off'. Default 'on'.
%
% 'displayInterval'     Number of iterations between plots.
%                       Default is 1 (plot after every iteration).
%
%======================================================================
% --Controlling reduction of dimension and whitening:
%
% Reduction of dimension is controlled by 'firstEig' and 'lastEig', or
% alternatively by 'interactivePCA'. 
%
% 'firstEig'            (integer) This and 'lastEig' specify the range for
%                       eigenvalues that are retained, 'firstEig' is
%                       the index of largest eigenvalue to be
%                       retained. Default is 1.
%
% 'lastEig'             (integer) This is the index of the last (smallest)
%                       eigenvalue to be retained. Default equals the
%                       dimension of data.
%
% 'interactivePCA'      (string) Either 'on' or 'off'. When set 'on', the
%                       eigenvalues are shown to the user and the
%                       range can be specified interactively. Default
%                       is 'off'.
%
% If you already know the eigenvalue decomposition of the covariance
% matrix, you can avoid computing it again by giving it with the
% following options:
%
% 'pcaE'                (matrix) Eigenvectors
% 'pcaD'                (matrix) Eigenvalues
%
% If you already know the whitened data, you can give it directly to
% the algorithm using the following options:
%
% 'whiteSig'            (matrix) Whitened signal
% 'whiteMat'            (matrix) Whitening matrix
% 'dewhiteMat'          (matrix) dewhitening matrix
%
% If values for all the 'whiteSig', 'whiteSig' and 'dewhiteMat' are
% suplied, they will be used in computing the ICA. PCA and whitening
% are not performed. Though 'mixedsig' is not used in the main
% algorithm it still must be entered - some values are still
% calculated from it.
%
% Performing preprocessing only is possible by the option:
%
% 'only'                (string) Compute only PCA i.e. reduction of
%                       dimension ('pca') or only PCA plus whitening
%                       ('white'). Default is 'all': do ICA estimation
%                       as well.  This option changes the output
%                       format accordingly. For example: 
%
%                       [whitesig, WM, DWM] = FASTICA(mixedsig, 
%                       'only', 'white') 
%                       returns the whitened signals, the whitening matrix
%                       (WM) and the dewhitening matrix (DWM). (See also
%                       WHITENV.) In FastICA the whitening matrix performs
%                       whitening and the reduction of dimension. Dewhitening
%                       matrix is the pseudoinverse of whitening matrix.
%                        
%                       [E, D] = FASTICA(mixedsig, 'only', 'pca') 
%                       returns the eigenvector (E) and diagonal 
%                       eigenvalue (D) matrices  containing the 
%                       selected subspaces. 
%
%======================================================================
% EXAMPLES
%
%       [icasig] = FASTICA (mixedsig, 'approach', 'symm', 'g', 'tanh');
%               Do ICA with tanh nonlinearity and in parallel (like
%               maximum likelihood estimation for supergaussian data).
%
%       [icasig] = FASTICA (mixedsig, 'lastEig', 10, 'numOfIC', 3);
%               Reduce dimension to 10, and estimate only 3
%               independent components.
%
%       [icasig] = FASTICA (mixedsig, 'verbose', 'off', 'displayMode', 'off');
%               Don't output convergence reports and don't plot
%               independent components.
%
%
% A graphical user interface for FASTICA can be launched by the
% command FASTICAG
%
%   See also FASTICAG

% 6.4.1998
% Hugo Gävert


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the mean and check the data

mixedsig = remmean(mixedsig);

[Dim, NumOfSampl] = size(mixedsig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values for optional parameters

% All
verbose           = 'on';

% Default values for 'pcamat' parameters
firstEig          = 1;
lastEig           = Dim;
interactivePCA    = 'off';

% Default values for 'fpica' parameters
approach          = 'defl';
numOfIC           = Dim;
g                 = 'pow3';
a1                = 1;
a2                = 1;
epsilon           = 0.0001;
maxNumIterations  = 1000;
initState         = 'rand';
guess             = 0;
displayMode       = 'on';
displayInterval   = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for fastICA - i.e. this file

b_verbose = 1;
jumpPCA = 0;
jumpWhitening = 0;
only = 3;
userNumOfIC = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the optional parameters

if(rem(nargin-1,2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:(nargin-1)/2
    % get the name and value of parameter
    str_param = eval (['s' int2str(i)]);
    val_param = eval (['v' int2str(i)]);

    % change the value of parameter
    if strcmp (str_param, 'verbose')
      verbose = val_param;
      % silence this program also
      if strcmp (verbose, 'off'), b_verbose = 0; end
    elseif strcmp (str_param, 'firstEig')
      firstEig = val_param;
    elseif strcmp (str_param, 'lastEig')
      lastEig = val_param;
    elseif strcmp (str_param, 'interactivePCA')
      interactivePCA = val_param;
    elseif strcmp (str_param, 'approach')
      approach = val_param;
    elseif strcmp (str_param, 'numOfIC')
      numOfIC = val_param;
      % User has suplied new value for numOfIC.
      % We'll use this information later on...
      userNumOfIC = 1;
    elseif strcmp (str_param, 'g')
      g = val_param;
    elseif strcmp (str_param, 'a1')
      a1 = val_param;
    elseif strcmp (str_param, 'a2')
      a2 = val_param;
    elseif strcmp (str_param, 'epsilon')
      epsilon = val_param;
    elseif strcmp (str_param, 'maxNumIterations')
      maxNumIterations = val_param;
    elseif strcmp (str_param, 'initGuess')
      % no use setting 'guess' if the 'initState' is not set
      initState = 'guess';
      guess = val_param;
    elseif strcmp (str_param, 'displayMode')
      displayMode = val_param;
    elseif strcmp (str_param, 'displayInterval')
      displayInterval = val_param;
    elseif strcmp (str_param, 'pcaE')
      % calculate if there are enought parameters to skip PCA
      jumpPCA = jumpPCA + 1;
      E = val_param;
    elseif strcmp (str_param, 'pcaD')
      % calculate if there are enought parameters to skip PCA
      jumpPCA = jumpPCA + 1;
      D = val_param;
    elseif strcmp (str_param, 'whiteSig')
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      whitesig = val_param;
    elseif strcmp (str_param, 'whiteMat')
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      whiteningMatrix = val_param;
    elseif strcmp (str_param, 'dewhiteMat')
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      dewhiteningMatrix = val_param;
    elseif strcmp (str_param, 'only')
      % if the user only wants to calculate PCA or...
      if strcmp(val_param, 'pca')
        only = 1;
      elseif strcmp(val_param, 'white')
        only = 2;
      elseif strcmp(val_param, 'all')
        only = 3;
      end

    else
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized parameter: ''' str_param '''']);
    end;
  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the data has been entered the wrong way,
% but warn only... it may be on purpose

if Dim > NumOfSampl
  if b_verbose
    fprintf('Warning: ');
    fprintf('The signal matrix may be oriented in the wrong way.\n');
    fprintf('In that case transpose the matrix.\n\n');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating PCA

% We need the results of PCA for whitening, but if we don't
% need to do whitening... then we dont need PCA...
if jumpWhitening == 3
  if b_verbose,
      fprintf ('Whitened signal and corresponding matrises suplied.\n');
      fprintf ('PCA calculations not needed.\n');
    end;
else

  % OK, so first we need to calculate PCA
  % Check to see if we already have the PCA data
  if jumpPCA == 2,
    if b_verbose,
      fprintf ('Values for PCA calculations suplied.\n');
      fprintf ('PCA calculations not needed.\n');
    end;
  else
    % display notice if the user entered one, but not both, of E and D.
    if (jumpPCA > 0) & (b_verbose),
      fprintf ('You must suply all of these in order to jump PCA:\n');
      fprintf ('''pcaE'', ''pcaD''.\n');
    end;

    % Calculate PCA
    [E, D]=pcamat(mixedsig, firstEig, lastEig, interactivePCA, verbose);
  end
end

% skip the rest if user only wanted PCA
if only > 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whitening the data

% Check to see if the whitening is needed...
if jumpWhitening == 3,
  if b_verbose,
    fprintf ('Whitening not needed.\n');
  end;
else

  % Whitening is needed
  % display notice if the user entered some of the whitening info, but not all.
  if (jumpWhitening > 0) & (b_verbose),
    fprintf ('You must suply all of these in order to jump whitening:\n');
    fprintf ('''whiteSig'', ''whiteMat'', ''dewhiteMat''.\n');
  end;

  % Calculate the whitening
  [whitesig, whiteningMatrix, dewhiteningMatrix] = ...
           whitenv (mixedsig, E, D, verbose);
end

end % if only > 1

% skip the rest if user only wanted PCA and whitening
if only > 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the ICA

% Check some parameters
% The dimension of the data may have been reduced during PCA calculations.
% The original dimension is calculated from the data by default, and the
% number of IC is by default set to equal that dimension.

Dim = size(whitesig, 1);

% In symmetric approach the number of IC must equal the dimension of data.
if strcmp (approach, 'symm')
  if numOfIC ~= Dim
    numOfIC = Dim;
    % Show warning only if verbose = 'on' and user suplied a value for 'numOfIC'
    if (b_verbose & userNumOfIC)
      fprintf('Warning: estimating %d independent components\n', numOfIC);
      fprintf('(In symmetric approach we estimate all available independent components)\n');
    end
  end
else
  % In deflation approach the number of IC can be different from the
  % dimension of data, but it must be less or equal.
  if numOfIC > Dim
    numOfIC = Dim;
    % Show warning only if verbose = 'on' and user suplied a value for 'numOfIC'
    if (b_verbose & userNumOfIC)
      fprintf('Warning: estimating only %d independent components\n', numOfIC);
      fprintf('(Can''t estimate more independent components than dimension of data)\n');
    end
  end
end

% Calculate the ICA with fixed point algorithm.
[A, W] = fpica (whitesig,  whiteningMatrix, dewhiteningMatrix, ...
  approach, numOfIC, g, a1, a2, epsilon, maxNumIterations, ...
  initState, guess, displayMode, displayInterval, verbose);


end % if only > 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output depends on the number of output parameters
% and the 'only' parameter.

if only == 1    % only PCA
  Out1 = E;
  Out2 = D;
elseif only == 2  % only PCA & whitening
  if nargout == 2
    Out1 = whiteningMatrix;
    Out2 = dewhiteningMatrix;
  else
    Out1 = whitesig;
    Out2 = whiteningMatrix;
    Out3 = dewhiteningMatrix;
  end
else      % ICA
  if nargout == 2
    Out1 = A;
    Out2 = W;
  else
    Out1 = W * mixedsig;
    Out2 = A;
    Out3 = W;
  end
end
