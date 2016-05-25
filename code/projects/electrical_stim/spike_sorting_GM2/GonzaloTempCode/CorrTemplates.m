function SimTem=CorrTemplates(templates,els)

for i=1:length(templates);
    for j=1:length(templates)
        SimTem(i,j)=trace(templates{i}(els,:)'*templates{j}(els,:));
    end
end