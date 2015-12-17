function plot_cell_bdry(B_use,cols,xlim1,ylim1)
 hold on;
 ncell=length(B_use);
    for icell=1:ncell
        
        ii=1:length(B_use{icell});
        %ii = convhull(B_use{icell}(:,1),B_use{icell}(:,2),'simplify',true);
    plot((B_use{icell}(ii,2)-ylim1(1)+1),(B_use{icell}(ii,1)-xlim1(1)+1),'Color',cols(icell+3,:),'LineWidth',2);
    end
    
end