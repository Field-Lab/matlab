   
function plot_SU_clustered(totalMaskAccept,label,su_log,nSU,total_mask_log)
[col_label]=distinguishable_colors(nSU+size(total_mask_log,2)+1+5);
        
[r,c] = find(totalMaskAccept);

x_coord = [min(r):max(r)];
y_coord  = [min(c):max(c)];

        mask = totalMaskAccept;
        sta_dim1 = size(mask,1);
        sta_dim2 = size(mask,2);
        indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
        masked_frame = indexedframe(logical(mask));
        
        % make plot
        ilist= 1:length(masked_frame);
        for ix=x_coord
            for iy=y_coord
                
                i1 =ilist(indexedframe(ix,iy)==masked_frame);
               
                if(~isempty(i1))
                  
                    for d1=[-1,0,1]
                        for d2=[-1,0,1]
                            if(d1==0 && d2==0)
                            continue;
                            end
                    i2 = ilist(indexedframe(ix+d1,iy+d2)==masked_frame);
                    
                    if(~isempty(i2))
                        
                            if(label(i1)==label(i2))
                                col=col_label(label(i1),:);
                            else
                                col=col_label(end,:);
                            end

                        if(su_log(i1,i2)>0)
                            plot([ix,ix+d1],[iy,iy+d2],'LineWidth',su_log(i1,i2)/4,'Color',col);
                            hold on;
%                             text(ix+d1/2,iy+d2/2,sprintf('%d',su_log(i1,i2)));
%                             hold on;
                        end
                    end
                    
                        end
                    end
                end
                
            end
        end
        
        for ix=x_coord
            for iy=y_coord
              i1 =ilist(indexedframe(ix,iy)==masked_frame);
                if(~isempty(i1))
                plot(ix,iy,'*','Color',col_label(label(i1),:));
                end
            end
        end
        
    axis equal   
    xlim([min(x_coord)-1,max(x_coord)+1]);
    ylim([min(y_coord)-1,max(y_coord)+1]);
ncell = size(total_mask_log,2);
B_use = cell(ncell,1);
B_use_ch = cell(ncell,1);
for icell=1:ncell
B = bwboundaries(reshape(total_mask_log(:,icell),[40,40]));
B_use{icell} = B{1};
B_use_ch{icell} = convhull(B{1}(:,1),B{1}(:,2),'simplify',true);
end


  hold on;
    for icell=1:ncell
        icell
        hold on;
    plot(B_use{icell}(B_use_ch{icell},1),B_use{icell}(B_use_ch{icell},2),'--','Color',col_label(icell+nSU,:),'LineWidth',3);
    end
   title(sprintf('# SU: %d',nSU));
   
end