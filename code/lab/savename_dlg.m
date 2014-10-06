function savename_dlg(handle)
% SAVENAME_DLG    Open a dialog window to set figure savename
%
% phli 2010-04
%

if nargin < 1
    handle = gcf;
end

answer = inputdlg('Savename:', 'Set savename', 1, {get_savename(handle)});

if length(answer) > 0
    set_savename(answer{1}, handle);
end
