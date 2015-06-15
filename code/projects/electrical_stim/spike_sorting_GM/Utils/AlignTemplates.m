function templates = AlignTemplates(templates)
%templates transform the  templates in a way so that all of them are aligned with onset of spike at time t=11.
%          It is assumed that for each neuron there is one that records spike at the soma,
%          which corresponds to the one that atains the minimum in the voltage values. That is the 
%          electrode that is aligned at t=11, and the rest of the electrodes are shifted accordingly
%input:     -templates, a cell array such that templates{n} is the template for neuron 
%          with index n, a matrix with E rows (one for each recording electrode).
%output     -templates: same templates but now with the
nNeurons= length(templates);
E = size(templates{1},1);
for n=1:nNeurons
    for e=1:E;
        temp=find(templates{n}(e,:)==min(templates{n}(e,:)));
        temp=temp(1);
        mins(e)=temp;
        minVal(e)=templates{n}(e,mins(e));
    end
    
    electrodeMin = find(minVal==min(minVal));
    electrodeMin = electrodeMin(1);
    
    dif = mins(electrodeMin)-11;
    
    if(dif>=0)
        Tmin = 1+dif;
        templatesNew{n}=templates{n}(:,Tmin:end);
    else

        templatesNew{n}(:,1:-dif)=repmat(templates{n}(:,1),1,-dif);
        templatesNew{n}(:,-dif+1:-dif+size(templates{n},2))=templates{n};
    end
end
templates=templatesNew;
