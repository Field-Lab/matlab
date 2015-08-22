function cone_peaks = inspect_cone_finding(all_sta, all_marks, final_sum, cone_peaks, base_diam, new_fit, coord_tform)

field_size = size(final_sum);


oneStep = 40;

 
for i=1:oneStep:field_size(1)
    x = i:min(field_size(1),i+oneStep);
    
    for j=1:oneStep:field_size(2)     
        y = j:min(field_size(2),j+oneStep);
        
        % find if there are cones in this region
        ifcones = find(cone_peaks(:,1)>=i & cone_peaks(:,1)<(i+oneStep) & cone_peaks(:,2)>=j & cone_peaks(:,2)<(j+oneStep));        
        if ~isempty(ifcones) 
            
            my_area = final_sum(x,y);
            
            figure(1)
            set(gcf,'position', [-1853         327        1027         774]);
            hold off
            colormap gray
            imagesc(my_area)
            hold on
            cone_handle = plot(cone_peaks(:,2)-j+1,cone_peaks(:,1)-i+1, 'xr');
            
            % find cells with marks in this area
            tmp = sum(squeeze(sum(all_marks(x,y,:), 1)));
%             tmp = sum(reshape(all_marks(x,y,:), length(x)*length(y), []));
            present_cells = find(tmp);
            
%             figure
%             imagesc(all_marks(:,:,present_cells(k)))
            
            
            [nrow,ncol] = opt_subplots(nnz(tmp));
            figure(2)
            set(gcf,'position', [-821   361   802   742]);
            ax = zeros(1,nnz(tmp));
            cone_handle_fig2_yel = ax;
            cone_handle_fig2_red = ax;
            for k = 1:nnz(tmp)
                ax(k) = subplot(nrow,ncol, k);
                hold off
                colormap gray
                imagesc(all_sta(:,:,present_cells(k)));
                hold on
                cone_handle_fig2_yel(k) = plot(cone_peaks(:, 2), cone_peaks(:, 1), 'y*');
                cone_handle_fig2_red(k) = plot(cone_peaks(ifcones, 2), cone_peaks(ifcones, 1), 'rx');
                
                [X, Y] = drawEllipse([new_fit.ctr(present_cells(k),:)*2 new_fit.rad(present_cells(k),:)*2 new_fit.fit_angle(present_cells(k))]);
                [X, Y] = tformfwd(coord_tform, X, Y);
                plot(X,Y,'r')
                line([j,j],[1,600], 'color', 'r')
                line([j,j]+oneStep,[1,600], 'color', 'r')
                line([1,600], [i,i],'color', 'r')
                line([1,600], [i,i]+oneStep,'color', 'r')
                
                axis([ y(1)-5, y(end)+5, x(1)-5,x(end)+5])
                set(gca,'dataaspectratio',[1 1 1])
                
            end
            linkaxes(ax)  
            
            but = 1;
            figure(1)
%             h = text(5,-1.5, 'READY', 'color','r', 'fontsize',24);
            while but == 1 || but == 3
                figure(1)
                pause(0.3)
                [click_xy(1), click_xy(2), but] = ginput(1);
                
%                 delete(h);
%                 h = text(5,-1.5, ['TOOK NEW!!  key = ', int2str(but)], 'color','r', 'fontsize',24);
                
                if but == 1 % if right click remove cone
                    
                    % find the closest                    
                    tmp=pdist2([click_xy(1),click_xy(2)],[cone_peaks(ifcones, 2)-j+1 cone_peaks(ifcones, 1)-i+1]);
                    [~,ind]=min(tmp);               
                    cone_peaks(ifcones(ind),:) = [];                    
                    ifcones(ind) = [];
                    ifcones(ind:end) = ifcones(ind:end)-1;
                elseif but==3 % if left click add cone
                    cone_peaks = [cone_peaks' [click_xy(2)+i-1; click_xy(1)+j-1]]';
                    ifcones =  [ifcones; size(cone_peaks,1)];
                end
                
                figure(1)
                delete(cone_handle);
                cone_handle = plot(cone_peaks(:,2)-j+1,cone_peaks(:,1)-i+1, 'xr');
                
                figure(2)
                for k=length(present_cells):-1:1
                    subplot(nrow,ncol,k)                    
                    delete(cone_handle_fig2_yel(k));
                    cone_handle_fig2_yel(k) = plot(cone_peaks(:, 2), cone_peaks(:, 1), 'y*');
                    if ishandle(cone_handle_fig2_red(k))
                        delete(cone_handle_fig2_red(k));
                    end
                    if ~isempty(ifcones)
                        cone_handle_fig2_red(k) = plot(cone_peaks(ifcones, 2), cone_peaks(ifcones, 1), 'rx');
                    end
                    disp('OK')
                end
                figure(1)
            end
            
            close(figure(2))
            
        end        
    end
end
% 
% for i=1:length(cone_peaks)
%     [X, Y] = drawEllipse([cone_peaks(i,:) [diam diam] pi/4]);
%     [X, Y] = tformfwd(coord_tform, X, Y);
% %     [inPoints] = polygrid( round(X), round(Y), 1);
% %     plot(inPoints(:,2), inPoints(:,1), 'x')
%     plot(round(Y),round(X),'r')
% end
% 


