function show_cones(flag, cellInd, varargin)

% index_type:   myCells - index of the cell in myCells array (from 1 to 5-10)
%               cellID - true index of cell (could be 3586)
%               datarun_id - index of cell in datarun (1 to N total cells)
%               cell_type - ord. index of cell in a cell_type array

p = inputParser;
p.addRequired('flag', @isnumeric);
p.addRequired('cellInd', @isnumeric);
p.addParamValue('index_type', 'myCells'); 
p.addParamValue('killCone', 1); 
% resolve user input and default values
p.parse(flag, cellInd, varargin{:});
% get params struct
params = p.Results;


global datarun myCells myLimits cones
persistent hConeInfo


% identify cell index in datarun
switch params.index_type
    case 'myCells'
        datInd=find(datarun.cell_ids==myCells(cellInd));
    case 'cellID'
        datInd=find(datarun.cell_ids==cellInd);
    case 'datarun_id'
        datInd=cellInd;
    case cell_type
        datInd=cell_indices(cellInd);
end

% get cell index in myCells array (for cones)
myInd=find(myCells==datarun.cell_ids(datInd));

% find subplot
hPlot=findobj('position',[0.4 0.62 0.3 0.3]);
subplot(hPlot)
hold on


if flag==1 % plot cones    

    for k=1:length(cones)
        if k~=myInd
            plot(cones{k}(:,1),cones{k}(:,2),'y+') % other cells
        else
            children = get(gca, 'children');
            for m=1:length(cones{k})
                myPoint=findobj('XData',cones{k}(m,1),'YData',cones{k}(m,2), 'color','r');
                if isempty(myPoint)
                    plot(cones{k}(m,1),cones{k}(m,2),'rx') % current cell - loop to make each datapoint deletable
                end
            end
        end
    end

elseif flag==2 % delete cone
    
    children = get(gca, 'children');
    myPoint=findobj('XData',cones{myInd}(params.killCone,1),'YData',cones{myInd}(params.killCone,2), 'color','r');
    delete(children(children==myPoint));
    cones{myInd}(params.killCone,:)=[];
        
end
    
%calculate new axis limits and adjust
new_x_min = min(myLimits(1), min(cones{myInd}(:,1))-5);
new_x_max = max(myLimits(2),max(cones{myInd}(:,1))+5);
new_y_min = min(myLimits(3),min(cones{myInd}(:,2))-5);
new_y_max = max(myLimits(4),max(cones{myInd}(:,2))+5);

axis([new_x_min new_x_max new_y_min new_y_max]);

% find largest ndd of cones
if size(cones{myInd},1)>1
    
    all_cones=cones{myInd};
    all_cones=squareform(pdist([all_cones(:,1),all_cones(:,2)]));
    all_cones(all_cones==0)=10000;
    all_cones=min(all_cones);
    mean_dist=mean(all_cones);
    std_dist=std(all_cones);
    cell_stat=['Mean nnd: ', num2str(mean_dist,2) ,' +/- ', num2str(std_dist,2)];
end


% print cell info
if ishandle(hConeInfo)
    delete(hConeInfo)
end

hConeInfo=uicontrol('style','text', 'Units', 'Normalized','position',[0.75 0.8 0.2 0.15],...
    'string',{['Cell ',int2str(myCells(myInd)), ': ', int2str(size(cones{myInd},1)),' cones'],'',...
    ['Largest nnd: ', num2str(max(all_cones),1)],'', cell_stat ,'',},...
    'fontsize',16, 'fontweight', 'bold');


% update overall cone info
update_cone_info;

show_all_cones;
