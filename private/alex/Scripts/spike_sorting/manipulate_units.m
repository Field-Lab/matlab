function manipulate_units(flag,xName,yName)
global rgc uicontr_units spikes2delete visibleUnits fatPoints signUnit ifUpdateISI ...
    subUnit selectValidUI

if isempty(rgc{1})
    flag=0;    
end

checked_units=zeros(1,size(rgc,2));
for i=1:size(rgc,2)
    checked_units(i)=get(uicontr_units(i),'Value');
end

switch flag
    case 1
        delete_unit;
        ifUpdateISI=1;
    case 2
        exclusive_spikes;
        ifUpdateISI=1;
    case 3
        join_unit;
        ifUpdateISI=1;
    case 4
        delete_spikes;
        ifUpdateISI=1;
    case 5
        visibleUnits=~checked_units;
    case 6
        fatPoints=checked_units;
    case 7
        signUnit=checked_units;
    case 8
        unit_from_overlap;
    case 9
        get_template(find(checked_units==1));
    case 10
        if get(selectValidUI,'Value')==1        
            if sum(checked_units)==0
                display('Check at least 1 unit')
            else
                visibleUnits=checked_units;
                if sum(checked_units)>1
                    subUnit=[rgc{logical(checked_units)}];
                    subUnit=unique(subUnit);
                else
                    subUnit=rgc{checked_units==1};
                end
            end
        else
            visibleUnits=ones(1,size(rgc,2));
        end
    case 11
        subset_pca(find(checked_units==1));
    otherwise
        display('There are no units so far')
        return
end



    function delete_unit
        keep_rgc=find(checked_units==0);
            if ~isempty(keep_rgc) % keep something
                tmp=1:size(rgc,2);
                tmp(keep_rgc)=[];
                rgc=rgc(keep_rgc);                
                visibleUnits(tmp)=[];
            else % delete everything
                rgc=cell(1,1);
                visibleUnits=[];
            end
    end

    function exclusive_spikes
        deleteWholeUnit=[];
        for j=find(checked_units==1)
%         for j=1:size(rgc,2)
%             if checked_units(j)
                for t=[1:j-1, j+1:length(rgc)]
                    [~,~,intersect_unit]=intersect(rgc{j},rgc{t});
                    if ~isempty(intersect_unit) % there is intersection
                        rgc{t}(intersect_unit)=[]; % remove spikes from another rgc    
                        if isempty(rgc{t}) % if all spikes deleted
                            deleteWholeUnit=[deleteWholeUnit t];
                        end
                    end
                end
%             end
        end
        deleteWholeUnit=unique(deleteWholeUnit);
        if ~isempty(deleteWholeUnit)
            tmp=1:size(rgc,2);
            tmp(deleteWholeUnit)=[];
            rgc=rgc(tmp);
            visibleUnits(deleteWholeUnit)=[];
        end
    end

    function join_unit
        joinUnits=find(checked_units==1);
        keepUnits=find(checked_units==0);
        if length(joinUnits)>1
            rgc{find(checked_units==1,1)}=[rgc{logical(checked_units)}];
            rgc{find(checked_units==1,1)}=unique(rgc{find(checked_units==1,1)});
%             for j=joinUnits(2:end)
%                 rgc{joinUnits(1)}=sort([rgc{joinUnits(1)}; rgc{j}]);
%             end
%             rgc{joinUnits(1)}=unique(rgc{joinUnits(1)});
            rgc=rgc(sort([keepUnits,joinUnits(1)]));
        end
    end

    function delete_spikes
        % units to delete spikes from
        deleteFrom=find(checked_units==1);
        toDelete=str2num(get(spikes2delete,'String')); % unit whose spikes to delete from others
        deleteWholeUnit=[];
        for j=deleteFrom
            [~,intersect_unit,~]=intersect(rgc{j},rgc{toDelete});
            if ~isempty(intersect_unit) % there is intersection
                rgc{j}(intersect_unit)=[]; % remove spikes from the rgc
            end
            if isempty(rgc{j}) % if all spikes deleted
                deleteWholeUnit=[deleteWholeUnit j];                
            end
        end
        if ~isempty(deleteWholeUnit)
            tmp=1:size(rgc,2);
            tmp(deleteWholeUnit)=[];
            rgc=rgc(tmp);
            visibleUnits(deleteWholeUnit)=[];
        end
    end

    function unit_from_overlap
        get_overlap_units=find(checked_units==1);
        if length(get_overlap_units)==2
            rgc{size(rgc,2)+1}=intersect(rgc{get_overlap_units(1)},rgc{get_overlap_units(2)});
            if isempty(rgc{end})
                display('No overlap!')
                rgc{end}=[];
            else
                visibleUnits(size(rgc,2))=1;
                fatPoints(size(rgc,2))=0;
                signUnit(size(rgc,2))=0;
                ifUpdateISI=1;
            end
        end
    end


redraw(xName,yName);
end