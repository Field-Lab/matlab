function set_savename(savename, fhandle)
% SET_SAVENAME    Set the UserData.savename for figure (defaults to gcf)
%
% usage: set_savename(savename, [fhandle])
%
% phli 2010-03
%

if nargin < 2
    fhandle = gcf;
end

set_ud(fhandle, 'savename', savename);