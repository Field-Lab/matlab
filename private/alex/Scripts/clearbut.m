function clearbut(allthere,leave)



for i=1:size(allthere,1)
    if strcmp(allthere(i).name,leave)
        eval(['clear ', allthere(i).name])
    end
end