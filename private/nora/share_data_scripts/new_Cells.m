count=0;
for i=1:56
    if sum(cells==OFF_Par(i))==0
        count=count+1;
        cells_new{count}=OFF_Par(i);
    else
        OFF_Par(i)
    end
end
for i=1:18
    if sum(cells==ON_Par(i))==0
        count=count+1;
        cells_new{count}=OFF_Par(i);
    else
        ON_Par(i)
    end
end
