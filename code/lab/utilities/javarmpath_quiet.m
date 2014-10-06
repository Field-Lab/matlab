function javarmpath_quiet(varargin)
% JAVARMPATH_QUIET  Remove from Java path without complaining if not found
% usage: javarmpath_quiet(varargin)
%
% Just wraps javarmpath, disabling the warning for things that are not
% found.  As an additional benefit, allows passing all the paths to remove
% as a cell array of strings in the first argument, in addition to the
% standard syntax of passing a varargin of strings.
%
% See also: JAVARMPATH
%
% phli 2010-05
%

if nargin == 1 && iscell(varargin{1})
    rms = varargin{1};
else
    rms = varargin;
end

warnquery = warning('query', 'MATLAB:GENERAL:JAVARMPATH:NotFoundInPath');
warning('off', 'MATLAB:GENERAL:JAVARMPATH:NotFoundInPath');
javarmpath(rms{:});
warning(warnquery.state, warnquery.identifier);