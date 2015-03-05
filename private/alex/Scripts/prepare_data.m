function [selected_data, cellIDs]=prepare_data(datarun, cell_list, nd, contr, condition, threshRev, threshStable, baselining, normalizing, ifReversable, ifStable)
% Parameters;
% datarun - data
% cell_list - off cells list or ON cells list
% contr - 'black' and 'white'
% threshRev - value for reversability check
% threshStable - threshold for stability check


switch nd
    case 6
        stableID=2;        
    case 5
        stableID=3;
    case 4
        stableID=4;        
end

switch condition
    case 'control'
        cond=1;
    case 'apb'
        cond=2;
    case 'wash'
        cond=3;
end

reverseID=stableID-1;

% cell list selection

% apply stability check
stableCellsList=intersect(cell_list, find(sum(squeeze(ifStable(:,stableID,:))>threshStable,2)==3));

% apply reversability check
cellIDs=intersect(stableCellsList,find(ifReversable(:,reverseID)>threshRev));

selected_data=[];
for i=cellIDs'
    data=eval(['datarun(',int2str(i),').',contr,'.nd',int2str(nd),'(:,',int2str(cond),')']);
    if baselining
        data=data-mean(data(1:450));
    end
    if normalizing %use OFF response peak for both flashes
        if strcmp(contr,'black')
            data=data/max(data(500:1000));
        else
            data=data/max(data(2500:3000));
        end
    end
    selected_data=[selected_data data];
end

