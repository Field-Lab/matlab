oncounter = 1;
offcounter = 1;

for j=1:length(cell_idx)
    if (mod(j,2) == 0 && offcounter <= length(offcell_idx))
        cell_idx(j) = offcell_idx(offcounter);
        fprintf('setting element %d to the %d th OFF cell\n',j,offcounter); 
        offcounter = offcounter + 1
        
    else
        if (mod(j,2) == 1 || offcounter > length(offcell_idx))
            cell_idx(j) = oncell_idx(oncounter);
            fprintf('setting element %d to the %d th ON cell\n',j,oncounter); 
            oncounter = oncounter + 1
        end
    end
    
end