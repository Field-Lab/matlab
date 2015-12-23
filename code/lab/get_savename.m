function savename = get_savename(fhandle)
% GET_SAVENAME    Get the figure savename, either from UserData or simply from date + figure number
%
% usage: get_savename([fhandle])
%        (defaults to gcf)
%
%
% 2010-03 phli
%

if nargin < 1
    fhandle = gcf;
end

savename = get_ud(fhandle, 'savename');
if isempty(savename)
   savename = [date '_fig' num2str(fhandle)];
end