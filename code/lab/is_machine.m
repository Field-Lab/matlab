function res=is_machine(name)
% is_machine('direwolf')
%
%greschner

[t,tt]=unix('hostname');
res=~isempty(strfind(tt,name));
