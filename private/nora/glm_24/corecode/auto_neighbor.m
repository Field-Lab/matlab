function [ON_Neighbors , OFF_Neighbors ] = auto_neighbor (init_pars , datarun, neighbornumber)

% New auth_neighbor on 2013-1-26, make neighbornumber an input param 

center = datarun{1}.vision.sta_fits{init_pars.master_idx}.mean;
center_x = center(1);
center_y = center(2);

% on oncolor = green
% off color  = magenta
MSHome  = 30;
MSNeigh = 20;
MSother = 8;
figure; hold on;
if strcmp ( init_pars.celltype , 'OFF-Parasol')
    plot(center_x , center_y, 'k.', 'markersize', MSHome); 
    plot(center_x , center_y, 'm*', 'markersize', MSHome);  
elseif strcmp (init_pars.celltype, 'ON-Parasol')
    plot(center_x , center_y, 'k.', 'markersize', MSHome); 
    plot(center_x , center_y, 'g*', 'markersize', MSHome); 
end    

for type = 1:2

    ids       = datarun{1}.cell_types{ type }.cell_ids;
    masteridx = zeros(size(ids));
    for i= 1:length(ids) 
        masteridx(i) = find (datarun{1}.cell_ids == ids(i) ) ;
    end
    RF      = datarun{1}.vision.sta_fits ( masteridx );
    
    %%%%%%%%%%%%%%%%%%%%%
    cells_loc = struct ('ids' ,ids, 'masteridx', masteridx, 'RF', RF);
    newcenters = zeros( length(ids) , 2);
    for i = 1:length(ids)
        newcntr = RF{i}.mean;
        newcenters(i,1) = newcntr(1);
        newcenters(i,2) = newcntr(2);
    end
    
    if type == 1
        plot( newcenters(: ,1) , newcenters(:,2) , 'g.' , 'markersize', MSother);
    elseif type ==2
        plot( newcenters(: ,1) , newcenters(:,2) , 'r.' , 'markersize', MSother);
    end
        
   
    for i = 1 : length(ids)
        mean = RF{i}.mean;       
        euclideandistance  = sqrt ( (mean(1)-center_x)^2  + (mean(2)-center_y)^2 );
        cells_loc(i).distance = euclideandistance;
    end
    
    distances = [cells_loc(:).distance];
    [a,b] = sort(distances,'ascend');
    

    if type == 1
        ON_Parasols = struct ('ids' ,ids, 'masteridx', masteridx, 'RF', RF);
        if a(1) < .001
            ON_Neighbors = ids(b(2:neighbornumber+1));
            shift = 1 ;
        else
            ON_Neighbors = ids(b(1:neighbornumber));
            shift = 0;
        end
        
        newcenters = zeros( length(ON_Neighbors) , 2);
    	for i = (1+shift):(shift+length(ON_Neighbors))
            newcntr = RF{b(i)}.mean;
            newcenters(i-shift,1) = newcntr(1);
            newcenters(i-shift,2) = newcntr(2);
        end
        plot( newcenters(: ,1) , newcenters(:,2) , 'k.' , 'markersize', MSNeigh);
        plot( newcenters(: ,1) , newcenters(:,2) , 'g*' , 'markersize', MSNeigh);        
    end
    
    
    if type == 2
        OFF_Parasols = struct ('ids' ,ids, 'masteridx', masteridx, 'RF', RF);
        if a(1) < .001
            OFF_Neighbors = ids(b(2:neighbornumber + 1));
            shift         = 1;
        else
            OFF_Neighbors = ids(b(1:neighbornumber));
            shift         = 0;
        end
        
        newcenters = zeros( length(OFF_Neighbors) , 2);
    	for i = (1+shift):(shift+length(OFF_Neighbors))
            newcntr = RF{b(i)}.mean;
            newcenters(i-shift,1) = newcntr(1);
            newcenters(i-shift,2) = newcntr(2);
        end
        plot( newcenters(: ,1) , newcenters(:,2) , 'k.' , 'markersize', MSNeigh);
        plot( newcenters(:,1) , newcenters(:,2) , 'r*' , 'markersize', MSNeigh);
    end
    
    
    
end
display( ON_Neighbors);
display(OFF_Neighbors);
hold off;

end









