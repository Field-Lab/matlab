function lock(fhandle)
% LOCK    Unset handle visibility to prevent overwriting figures
%
% usage: lock([fhandle])
%
% inputs: fhandle - Optional figure handle; defaults to gcf
%
% example: lock % Sets gcf HandleVisibility to 'off'; figure no longer gets overwritten
%
% 2010-02 phli
%

if nargin < 1
    fhandle = gcf;
end

set(fhandle, 'HandleVisibility', 'off');