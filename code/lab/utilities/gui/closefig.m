function closefig(h)
% CLOSEFIG     Close the figure that contains the object referred to by handle H
%
% usage:  closefig(axes_or_figure_handle)
%
% arguments:  h - Handle for the object whose ancestor figure to close
%
% 2010-01 phli

close(getfig(h));