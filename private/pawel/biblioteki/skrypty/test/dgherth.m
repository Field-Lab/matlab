All=[];
for i=1:64
    n=PatternsForMovie(i);
    name=['Pattern_' num2str(n)]
    f=fopen('Pattern_397');
    a=fread(f);
    b=reshape(a,2,length(a)/2);
    el=b(1,:);
    electrodes=find(el>255);
    All=[All el]
    fclose(f);
end