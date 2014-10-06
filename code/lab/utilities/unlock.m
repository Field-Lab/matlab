function unlock(fhandle)
% UNLOCK    Set handle visibility to allow overwriting figures
%
% usage: unlock(fhandle)
%
% inputs: fhandle - Figure handle; must be given because if a figure is
%                   already locked it will not show up with gcf.
%
%
% 2010-02 phli
%

set(fhandle, 'HandleVisibility', 'on');