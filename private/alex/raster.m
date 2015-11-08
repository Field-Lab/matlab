
clear

NoF=739;
%flicks
file_types=zeros(1,NoF);
file_types([3+2:3:540])=1;
file_types([544+2:3:NoF])=1;

file_list=find(file_types(1:NoF)==1);


t=[];
cnt=0;
for i=file_list
    t=[t unit{i}*1000+62000*cnt];
    cnt=cnt+1;
end
rasterplot(t,cnt,62000)