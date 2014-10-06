function the_name = printable_name(datarun)
% printable_name     Return name of the datarun which can be printed in a matlab figure
%
% usage:  the_name = printable_name(datarun)
%
% arguments:     datarun - datarun struct
%
% outputs:     the_name - string
%
%
% 2008-10 gauthier
%

if isfield(datarun,'names') && isfield(datarun.names,'short_name') && ~isempty(datarun.names.short_name)
    the_name = strrep(datarun.names.short_name,'_',' / ');
else
    the_name = 'unknown datarun';
end


