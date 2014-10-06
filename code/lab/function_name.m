function function_name
%prints name of calling functions on top of plot
%greschner

t=dbstack;
name=[];
for i=length(t):-1:2
    name=[name '  ' t(i).name '(' num2str(t(i).line) ')'];
end
name=strrep(name,'_','\_');

%position=get(gca,'position')
%visible=get(gca,'visible');

axes('position',[0,0,1,1],'visible','off');

text(0.13,0.94,name,'FontSize',8);

%axes('position',position,'visible',visible);