function handle_cones(flag, cellInd, varargin)

% index_type:   myCells - index of the cell in myCells array (from 1 to 5-10)
%               cellID - true index of cell (could be 3586)
%               datarun_id - index of cell in datarun (1 to N total cells)
%               cell_type - ord. index of cell in a cell_type array

global datarun myCells cones STAplotPosition

p = inputParser;
p.addRequired('flag', @isnumeric);
p.addRequired('cellInd', @isnumeric);
p.addParamValue('index_type', 'myCells');
p.addParamValue('frame', 5);
p.addParamValue('nnd_scale', 3);
% resolve user input and default values
p.parse(flag, cellInd, varargin{:});
% get params struct
params = p.Results;


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
hPlot=findobj('position',STAplotPosition);
subplot(hPlot)

if flag==0 % find cones automatically
    
    % prepare sta
    sta=datarun.stas.stas{datInd};
    sta=sta(:,:,:,params.frame); % frame
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
    tmp_sta=sta;
    if abs(min(tmp_sta(:)))>max(tmp_sta(:))  % OFF cell, invert polarity
        tmp_sta=-tmp_sta;
    end
    
    
    if length(cones)<myInd % first time, find a bunch
        cones{myInd}=[]; 

        w_center=[];

        keep_looking=true;
        while keep_looking && length(cones{myInd})<7
            myMax=max(tmp_sta(:));
            myCoord=find(tmp_sta==myMax,1);
            [row, col]=ind2sub(size(tmp_sta),myCoord);
            % check if border position... if, then correct - don't trust!
            row=min(max(row,3),size(sta,1)-3);
            col=min(max(col,3),size(sta,2)-3);
            
            if length(cones{myInd})>5 % start checking from 6th cone, 5 are obligatory
                
                % center of mass: check if new cone is not at most 2
                % of max distance between center and other cones
                x_center=cones{myInd}(:,1);
                y_center=cones{myInd}(:,2);
                x_com=double(sum(w_center.*x_center)/sum(w_center));
                y_com=double(sum(w_center.*y_center)/sum(w_center));
                
                all_dist=pdist2([x_com y_com], [x_center,y_center]);
                max_dist=max(all_dist);
                
                if pdist2([col row], [x_com y_com])>2*max_dist
                    keep_looking=false;
                end
                
                % mean nnd
                tmp=squareform(pdist([x_center,y_center]));
                tmp(tmp==0)=10000;
                tmp=min(tmp);
                if min(pdist2([col row], [x_center,y_center]))>(mean(tmp)+params.nnd_scale*std(tmp))
                    keep_looking=false;
                end
            end
            if keep_looking
                w_center=[w_center; tmp_sta(row,col)];
                tmp_sta(row-2:row+2,col-2:col+2,:)=0;
                cones{myInd}=[cones{myInd}; [col, row]];
            end
        end
        
        delete_cell(0); % check if this cone is on deleted list; remove and correct display
        
    else % automatically add cones (1 by 1)
        
        for i=1:length(cones{myInd})
            col=cones{myInd}(i,1);
            row=cones{myInd}(i,2);
            tmp_sta(row-2:row+2,col-2:col+2,:)=0;
        end
        
        myMax=max(tmp_sta(:));
        myCoord=find(tmp_sta==myMax,1);
        [row, col]=ind2sub(size(tmp_sta),myCoord);
        % check if border position... if, then correct - don't trust!
        row=min(max(row,3),size(sta,1)-3);
        col=min(max(col,3),size(sta,2)-3);        
        cones{myInd}=[cones{myInd}; [col, row]];
        
    end
    
    show_cones(1,myInd,'index_type','myCells');
    
elseif flag==1 % add manually
    
    
    while 1
        
        subplot(hPlot);
        plot_position=get(hPlot,'position');
        
        Xlim=get(gca,'Xlim');
        Ylim=get(gca,'Ylim');
        
        [x,y]=ginput(1);
        clicked_plot_position=get(gca,'position');    
        dd=nnz(plot_position-clicked_plot_position);
       
        if dd || x<Xlim(1) || x>Xlim(2) || y<Ylim(1) || y>Ylim(2)
            break;
        end
        
        cones{myInd}=[cones{myInd}; [x, y]];
        
        show_cones(1,myInd,'index_type','myCells')
        
    end
    
elseif flag==2 % delete manually
    
    while 1
        
        subplot(hPlot);
        plot_position=get(hPlot,'position');
        
        Xlim=get(gca,'Xlim');
        Ylim=get(gca,'Ylim');
        
        [x,y]=ginput(1);
        clicked_plot_position=get(gca,'position');    
        dd=nnz(plot_position-clicked_plot_position);
       
        if dd || x<Xlim(1) || x>Xlim(2) || y<Ylim(1) || y>Ylim(2)
            break;
        end
        
        [~, ind]=min(pdist2([x y], cones{myInd})); 
        
        show_cones(2,myInd,'index_type','myCells','killCone',ind)

    end
    
elseif flag==3 % delete automatically, last cone in the current cell
    
    ind=length(cones{myInd});
    
    show_cones(2,myInd,'index_type','myCells','killCone',ind)
    
elseif flag==4 % delete all cones in this cell, but don't delete the cell
    
        show_cones(4,myInd,'index_type','myCells','killCone',length(cones{myInd}))
    
end

