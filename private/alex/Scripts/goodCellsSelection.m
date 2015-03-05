function [selected_data, cellIDs]=goodCellsSelection(datarun, cell_list, nd, contr, condition,varargin)
% goodCellsSelection     selects cells and prepares them for analysis
%
% usage:   [selected_data, cellIDs]=goodCellsSelection(datarun, cell_list, nd, contr, condition, <params>)
%
% arguments:    datarun - datarun struct
%               cell_list - list of cells preselected from datarun (OFF cells in particular)
%               nd - ND, numeric 6, 5, or 4
%               contr - contrast of the stimulus,string, 'black' or 'white'
%               condition - string,'control', 'apb', or 'wash'
%               params - struct or list of optional parameters (see below)
%
% outputs:      selected data - processed (quality checked, baselined, normalized) responses of selected cells
%               cellIDs - list of cells selected after quality checks
%
% parameters:   'threshRev', threshold for 'reversability check', -1 to 1, default:
%                   -1, will take all cells
%               'threshStable', threshold for 'stability check', -1 to 1, default:
%                   -1, will take all cells
%               'baselining', boolean, default: true
%               'normalizing', boolean, default: true



p = inputParser;

% specify list of optional parameters
p.addParamValue('threshRev', -1);
p.addParamValue('threshStable', -1);
p.addParamValue('baselining', true);
p.addParamValue('normalizing', true);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



if strcmp(contr,'black')
    contrastCorrect=0;
else
    contrastCorrect=2000;
end


switch condition
    case 'control'
        cond=1;
    case 'apb'
        cond=2;
    case 'wash'
        cond=3;
end

% reversible cells
if params.threshRev>-1
    ranges=500:2000;
    reverseList=[];
    for i=cell_list
        data1=eval(['datarun(i).',contr,'.nd',int2str(nd),'(ranges+contrastCorrect,1)']); % control
        data2=eval(['datarun(i).',contr,'.nd',int2str(nd),'(ranges+contrastCorrect,3)']); % wash
        if corr(data1,data2)>params.threshRev
               reverseList=[reverseList i];
        end
    end
else 
    reverseList=cell_list;
end

%stable cells
if params.threshStable>-1
    stableList=[];
    for i=cell_list
        data=eval(['datarun(i).consistency.', contr, '.nd',int2str(nd),'(',int2str(cond),')']);
        if data>params.threshStable
            stableList=[stableList i];
        end
    end
else
    stableList=cell_list;
end


% quality checks
cellIDs=intersect(cell_list, intersect(stableList,reverseList));


selected_data=[];
for i=cellIDs
    data=eval(['datarun(i).',contr,'.nd',int2str(nd),'(:,',int2str(cond),')']);
    if params.baselining
        data=data-mean(data(1:450));
    end
    if params.normalizing %use OFF response peak for both flashes
        if strcmp(contr,'black')
            data=data/max(data(500:1000));
        else
            data=data/max(data(2500:3000));
        end
    end
    selected_data=[selected_data data];
end

