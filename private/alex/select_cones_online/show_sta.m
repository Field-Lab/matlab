function show_sta(cellInd, varargin)

% index_type:   myCells - index of the cell in myCells array (from 1 to 5-10)
%               cellID - true index of cell (could be 3586)
%               datarun_id - index of cell in datarun (1 to N total cells)
%               cell_type - ord. index of cell in a cell_type array
% scale:        scale factor to show RF fit in units of max radius
% frame:        frame of sta - button to change

global datarun myCells myLimits cell_indices coord_tform ctr rad fit_angle
persistent hPlot

p = inputParser;
p.addRequired('cellInd', @isnumeric);
p.addParamValue('scale', 2); 
p.addParamValue('index_type', 'myCells'); 
p.addParamValue('frame', 5); 

% resolve user input and default values
p.parse(cellInd, varargin{:});

% get params struct
params = p.Results;

hPlot=findobj('position',[0.4 0.62 0.3 0.3]);

if cellInd>0 % if not 0, show sta

    if ~isempty(hPlot)
        delete(hPlot)
    end
    hPlot=subplot('position',[0.4 0.62 0.3 0.3]);
    set(hPlot,'DataAspectRatio',[1 1 1]);
    axis ij
    
    
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
    
    % check and load sta if not loaded
    if isempty(datarun.stas.stas{datInd})
        datarun = load_sta(datarun, struct('load_sta',myCells(myInd)));
    end
    
    % prepare sta
    sta=squeeze(datarun.stas.stas{datInd});
    sta=sta(:,:,params.frame); % frame
    
    % prepare RF fit
    centre=round(ctr(datInd,:));
    radius=ceil(max(rad(datInd,:)))*params.scale;
    [X, Y] = drawEllipse([ctr(datInd,:) rad(datInd,:) fit_angle(datInd)]);
    [X, Y] = tformfwd(coord_tform, X, Y);
    
    % plot stuff
    subplot(hPlot);
    colormap gray
    
    imagesc(sta);
    
    hold on
    
    plot(X,Y,'r')
    
    % get axis limits
    minx=max(centre(1)-radius,1);
    maxx=min(centre(1)+radius,datarun.stimulus.field_width);
    miny=max(centre(2)-radius,1);
    maxy=min(centre(2)+radius, datarun.stimulus.field_height);
    
    myLimits=[minx maxx miny maxy];

    axis(myLimits)
    
elseif cellInd==0 % zoom out by 30% each time
    
    subplot(hPlot)
    xx=get(hPlot,'XLim');
    yy=get(hPlot,'YLim');
    
    tmpx=diff(xx)*0.3;
    tmpy=diff(yy)*0.3;
    xx(1)=max(xx(1)-tmpx,1);    
    xx(2)=min(xx(2)+tmpx,datarun.stimulus.field_width);    
    yy(1)=max(yy(1)-tmpy,1);
    yy(2)=min(yy(2)+tmpy,datarun.stimulus.field_height);
    
    axis([xx yy])
    
    
end


