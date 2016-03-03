function sta = show_sta(cellInd, varargin)

% index_type:   myCells - index of the cell in myCells array (from 1 to 5-10)
%               cellID - true index of cell (could be 3586)
%               datarun_id - index of cell in datarun (1 to N total cells)
%               cell_type - ord. index of cell in a cell_type array
% scale:        scale factor to show RF fit in units of max radius
% frame:        frame of sta - button to change

global datarun myCells myLimits cell_indices coord_tform ctr rad fit_angle STAplotPosition
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

hPlot=findobj('position',STAplotPosition);

if cellInd>0 % if not 0, show sta

    if ~isempty(hPlot)
        delete(hPlot)
    end
    hPlot=subplot('position',STAplotPosition);
%     set(hPlot,'DataAspectRatio',[1 1 1]);
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
    
    % prepare sta
    sta=datarun.stas.stas{datInd};
    aa=sta(:,:,:,params.frame);
    bb=sta(:,:,:,params.frame-1);
    
    if max(abs(aa(:)))<max(abs(bb(:)));
        sta = bb;
    else
        sta=sta(:,:,:,params.frame); % frame
    end
%     sta=sta(:,:,:,params.frame); % frame
    if datarun.stimulus.independent=='t'
        % RGB run, choose either green or blue channel
        tmpGreen=sta(:,:,2);
        tmpBlue=sta(:,:,3);
        if max(abs(tmpGreen(:)))>max(tmpBlue(:)) % green gun bigger; non Blue cell
            sta=tmpGreen;
        else
            sta=tmpBlue;
        end
    end
    sta=double(squeeze(sta));
    
    subplot(hPlot);
    colormap gray
    
    imagesc(sta);
 
end


