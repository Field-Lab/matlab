function plot_cone_interaction_fits(params_ref_all, params_x_all, params_y_all, ...
    params_xy_all, resnorm_x_all, resnorm_y_all, resnorm_xy_all, x_all, y_all, ...
    center_cones, path2save, datarunID, cellID)


for cone1 = 1:length(center_cones)-1
    
    for cone2=cone1+1:length(center_cones)
        
        x2 = x_all{cone1};
        y2 = y_all{cone1, cone2};
        params_ref = params_ref_all{cone1, cone2};
        params_x = params_x_all{cone1, cone2};
        params_y = params_y_all{cone1, cone2};
        params_xy = params_xy_all{cone1, cone2};
        resnorm_x = resnorm_x_all(cone1,cone2)*100;
        resnorm_y = resnorm_y_all(cone1,cone2)*100;
        resnorm_xy = resnorm_xy_all(cone1,cone2)*100;
        
        
        x22 = x_all{cone2};
        y22 = y_all{cone2, cone1};
        params_ref2 = params_ref_all{cone2, cone1};
        params_x2 = params_x_all{cone2, cone1};
        params_y2 = params_y_all{cone2, cone1};
        params_xy2 = params_xy_all{cone2, cone1};
        resnorm_x2 = resnorm_x_all(cone2,cone1)*100;
        resnorm_y2 = resnorm_y_all(cone2,cone1)*100;
        resnorm_xy2 = resnorm_xy_all(cone2,cone1)*100;
        
        figure
        set(gcf, 'position', [-1703 260 981 789], 'Name', ['Ref cone ', int2str(cone1), ', cond cone ', int2str(cone2)])
        
        
        subplot(3,3,1)
        x = x2(:,1);
        sat   = params_ref(1);
        sigma = params_ref(2);
        mu = params_ref(3);
        sh = params_ref(4);
        y = sat .* normcdf(x, mu, sigma)+sh;
        hold on
        plot(x,y, '-*')
        plot(x2,y2)
        title(['original cone ', int2str(center_cones(cone1))])
        
        
        subplot(3,3,2)
        x = x22(:,1);
        sat   = params_ref2(1);
        sigma = params_ref2(2);
        mu = params_ref2(3);
        sh = params_ref2(4);
        y = sat .* normcdf(x, mu, sigma)+sh;
        hold on
        plot(x,y, '-*')
        plot(x22,y22)
        title(['original cone ', int2str(center_cones(cone2))])
        
        subplot(3,3,3)
        hold on
        plot(resnorm_x_all(cone1,:)*100, resnorm_y_all(cone1,:)*100, 'ob')
        plot(resnorm_x_all(cone2,:)*100, resnorm_y_all(cone2,:)*100, 'or')
        plot(resnorm_x_all(cone1,cone2)*100, resnorm_y_all(cone1,cone2)*100, 'xk')
        plot(resnorm_x_all(cone2,cone1)*100, resnorm_y_all(cone2,cone1)*100, '+k')
        tmp = max([resnorm_x_all(cone1,:)*100 resnorm_y_all(cone1,:)*100 ...
            resnorm_x_all(cone2,:)*100 resnorm_y_all(cone2,:)*100]);
        line([0,tmp], [0 tmp], 'color', 'k')
        axis tight
        xlabel('x shift')
        ylabel('y shift')
        title(['blue ', int2str(center_cones(cone1)) ', red ', int2str(center_cones(cone2))])
        
        subplot(3,3,4)
        x = x2(1,1)+min(params_x(1:3)):0.01:x2(end,1)+max(params_x(1:3));
        n = size(x2,2);
        sat   = params_x(n+1);
        sigma = params_x(n+2);
        mu = 0;
        sh = params_x(n+3);
        y = sat .* normcdf(x, mu, sigma)+sh;
        hold on
        plot(x,y, '-*')
        hold on
        plot(x2(:,1)+params_x(1),y2(:,1))
        plot(x2(:,1)+params_x(2),y2(:,2))
        plot(x2(:,1)+params_x(3),y2(:,3))
        title(['x shift, ', num2str(resnorm_x)])
        
        subplot(3,3,5)
        x = x2(:,1);
        n = size(x2,2);
        sat   = params_y(n+1);
        sigma = params_y(n+2);
        mu = params_y(n+3);
        y = sat .* normcdf(x, mu, sigma);
        hold on
        plot(x,y, '-*')
        plot(x,y2(:,1)-params_y(1))
        plot(x,y2(:,2)-params_y(2))
        plot(x,y2(:,3)-params_y(3))
        title(['y shift, ', num2str(resnorm_y)])
        
        subplot(3,3,6)
        x = x2(1,1)+min(params_xy(2:2:6)):0.01:x2(end,1)+max(params_xy(2:2:6));
        n = size(x2,2)*2;
        sat   = params_xy(n+1);
        sigma = params_xy(n+2);
        mu = 0;
        y = sat .* normcdf(x, mu, sigma);
        hold on
        plot(x,y, '-*')
        plot(x2(:,1)+params_xy(2),y2(:,1)-params_xy(1))
        plot(x2(:,1)+params_xy(4),y2(:,2)-params_xy(3))
        plot(x2(:,1)+params_xy(6),y2(:,3)-params_xy(5))
        title(['xy shift, ', num2str(resnorm_xy)])
        
        
        subplot(3,3,7)
        x = x22(1,1)+min(params_x2(1:3)):0.01:x22(end,1)+max(params_x2(1:3));
        n = size(x22,2);
        sat   = params_x2(n+1);
        sigma = params_x2(n+2);
        mu = 0;
        sh = params_x2(n+3);
        y = sat .* normcdf(x, mu, sigma)+sh;
        hold on
        plot(x,y, '-*')
        hold on
        plot(x22(:,1)+params_x2(1),y22(:,1))
        plot(x22(:,1)+params_x2(2),y22(:,2))
        plot(x22(:,1)+params_x2(3),y22(:,3))
        title(['x shift, ', num2str(resnorm_x2)])
        
        subplot(3,3,8)
        x = x22(:,1);
        n = size(x22,2);
        sat   = params_y2(n+1);
        sigma = params_y2(n+2);
        mu = params_y2(n+3);
        y = sat .* normcdf(x, mu, sigma);
        hold on
        plot(x,y, '-*')
        plot(x,y22(:,1)-params_y2(1))
        plot(x,y22(:,2)-params_y2(2))
        plot(x,y22(:,3)-params_y2(3))
        title(['y shift, ', num2str(resnorm_y2)])
        
        subplot(3,3,9)
        x = x22(1,1)+min(params_xy2(2:2:6)):0.01:x22(end,1)+max(params_xy2(2:2:6));
        n = size(x2,2)*2;
        sat   = params_xy2(n+1);
        sigma = params_xy2(n+2);
        mu = 0;
        y = sat .* normcdf(x, mu, sigma);
        hold on
        plot(x,y, '-*')
        plot(x22(:,1)+params_xy2(2),y22(:,1)-params_xy2(1))
        plot(x22(:,1)+params_xy2(4),y22(:,2)-params_xy2(3))
        plot(x22(:,1)+params_xy2(6),y22(:,3)-params_xy2(5))
        title(['xy shift, ', num2str(resnorm_xy2)])

        drawnow             
        if ~isdir([path2save,int2str(datarunID), '/tiff/'])
            mkdir([path2save,int2str(datarunID), '/tiff/'])
        end
        saveas(gcf,[path2save,int2str(datarunID), '/tiff/ID_',int2str(cellID),'_cones_', int2str(cone1),'_', int2str(cone2),'.tiff'])
        
        close(gcf)
    end
    
end


% if 0
% 
% for cone1 = 1:length(center_cones)-1
%     
%     for cone2=cone1+1:length(center_cones)
%         
%         x2 = x_all{cone1};
%         y2 = y_all{cone1, cone2};
%         params_ref = params_ref_all{cone1, cone2};
%         params_x = params_x_all{cone1, cone2};
%         params_y = params_y_all{cone1, cone2};
%         params_xy = params_xy_all{cone1, cone2};
%         resnorm_x = resnorm_x_all(cone1,cone2)*100;
%         resnorm_y = resnorm_y_all(cone1,cone2)*100;
%         resnorm_xy = resnorm_xy_all(cone1,cone2)*100;
%         
%         
%         x22 = x_all{cone2};
%         y22 = y_all{cone2, cone1};
%         params_ref2 = params_ref_all{cone2, cone1};
%         params_x2 = params_x_all{cone2, cone1};
%         params_y2 = params_y_all{cone2, cone1};
%         params_xy2 = params_xy_all{cone2, cone1};
%         resnorm_x2 = resnorm_x_all(cone2,cone1)*100;
%         resnorm_y2 = resnorm_y_all(cone2,cone1)*100;
%         resnorm_xy2 = resnorm_xy_all(cone2,cone1)*100;
%         
%         figure
%         set(gcf, 'position', [-1703 260 981 789], 'Name', ['Ref cone ', int2str(cone1), ', cond cone ', int2str(cone2)])
%         
%         
%         subplot(3,3,1)
%         x = x2(:,1);
%         sat   = params_ref(1);
%         sigma = params_ref(2);
%         mu = params_ref(3);
%         sh = params_ref(4);
%         y = sat .* normcdf(x, mu, sigma)+sh;
%         hold on
%         plot(x,y, '-*')
%         plot(x2,y2)
%         title(['original cone ', int2str(center_cones(cone1))])
%         
%         
%         subplot(3,3,2)
%         x = x22(:,1);
%         sat   = params_ref2(1);
%         sigma = params_ref2(2);
%         mu = params_ref2(3);
%         sh = params_ref2(4);
%         y = sat .* normcdf(x, mu, sigma)+sh;
%         hold on
%         plot(x,y, '-*')
%         plot(x22,y22)
%         title(['original cone ', int2str(center_cones(cone2))])
%         
%         subplot(3,3,3)
%         hold on
%         plot(resnorm_x_all(cone1,:)*100, resnorm_y_all(cone1,:)*100, 'ob')
%         plot(resnorm_x_all(cone2,:)*100, resnorm_y_all(cone2,:)*100, 'or')
%         plot(resnorm_x_all(cone1,cone2)*100, resnorm_y_all(cone1,cone2)*100, 'xk')
%         plot(resnorm_x_all(cone2,cone1)*100, resnorm_y_all(cone2,cone1)*100, '+k')
%         tmp = max([resnorm_x_all(cone1,:)*100 resnorm_y_all(cone1,:)*100 ...
%             resnorm_x_all(cone2,:)*100 resnorm_y_all(cone2,:)*100]);
%         line([0,tmp], [0 tmp], 'color', 'k')
%         axis tight
%         xlabel('x shift')
%         ylabel('y shift')
%         title(['blue ', int2str(center_cones(cone1)) ', red ', int2str(center_cones(cone2))])
%         
%         subplot(3,3,4)
%         x = x2(1,1)+min(params_x(1:3)):0.05:x2(end,1)+max(params_x(1:3));
%         n = size(x2,2);
%         sat   = params_x(n+1);
%         sigma = params_x(n+2);
%         mu = params_x(n+3);
%         sh = params_x(n+4);
%         y = sat .* normcdf(x, mu, sigma)+sh;
%         hold on
%         plot(x,y, '-*')
%         plot(x2(:,1)+params_x(1),y2(:,1))
%         plot(x2(:,1)+params_x(2),y2(:,2))
%         plot(x2(:,1)+params_x(3),y2(:,3))
%         title(['x shift, ', num2str(resnorm_x)])
%         
%         subplot(3,3,5)
%         x = x2(:,1);
%         n = size(x2,2);
%         sat   = params_y(n+1);
%         sigma = params_y(n+2);
%         mu = params_y(n+3);
%         sh = params_y(n+4);
%         y = sat .* normcdf(x, mu, sigma)+sh;
%         hold on
%         plot(x,y, '-*')
%         plot(x,y2(:,1)-params_y(1))
%         plot(x,y2(:,2)-params_y(2))
%         plot(x,y2(:,3)-params_y(3))
%         title(['y shift, ', num2str(resnorm_y)])
%         
%         subplot(3,3,6)
%         x = x2(1,1)+min(params_xy(2:2:6)):0.05:x2(end,1)+max(params_xy(2:2:6));
%         n = size(x2,2)*2;
%         sat   = params_xy(n+1);
%         sigma = params_xy(n+2);
%         mu = params_xy(n+3);
%         sh = params_xy(n+4);
%         y = sat .* normcdf(x, mu, sigma)+sh;
%         hold on
%         plot(x,y, '-*')
%         plot(x2(:,1)+params_xy(2),y2(:,1)-params_xy(1))
%         plot(x2(:,1)+params_xy(4),y2(:,2)-params_xy(3))
%         plot(x2(:,1)+params_xy(6),y2(:,3)-params_xy(5))
%         title(['xy shift, ', num2str(resnorm_xy)])
%         
%         
%         subplot(3,3,7)
%         x = x22(1,1)+min(params_x2(1:3)):0.05:x22(end,1)+max(params_x2(1:3));
%         n = size(x22,2);
%         sat   = params_x2(n+1);
%         sigma = params_x2(n+2);
%         mu = params_x2(n+3);
%         sh = params_x2(n+4);
%         y = sat .* normcdf(x, mu, sigma)+sh;
%         hold on
%         plot(x,y, '-*')
%         hold on
%         plot(x22(:,1)+params_x2(1),y22(:,1))
%         plot(x22(:,1)+params_x2(2),y22(:,2))
%         plot(x22(:,1)+params_x2(3),y22(:,3))
%         title(['x shift, ', num2str(resnorm_x2)])
%         
%         subplot(3,3,8)
%         x = x22(:,1);
%         n = size(x22,2);
%         sat   = params_y2(n+1);
%         sigma = params_y2(n+2);
%         mu = params_y2(n+3);
%         sh = params_y2(n+4);
%         y = sat .* normcdf(x, mu, sigma)+sh;
%         hold on
%         plot(x,y, '-*')
%         plot(x,y22(:,1)-params_y2(1))
%         plot(x,y22(:,2)-params_y2(2))
%         plot(x,y22(:,3)-params_y2(3))
%         title(['y shift, ', num2str(resnorm_y2)])
%         
%         subplot(3,3,9)
%         x = x22(1,1)+min(params_xy2(2:2:6)):0.05:x22(end,1)+max(params_xy2(2:2:6));
%         n = size(x2,2)*2;
%         sat   = params_xy2(n+1);
%         sigma = params_xy2(n+2);
%         mu = params_xy2(n+3);
%         sh = params_xy2(n+4);
%         y = sat .* normcdf(x, mu, sigma)+sh;
%         hold on
%         plot(x,y, '-*')
%         plot(x22(:,1)+params_xy2(2),y22(:,1)-params_xy2(1))
%         plot(x22(:,1)+params_xy2(4),y22(:,2)-params_xy2(3))
%         plot(x22(:,1)+params_xy2(6),y22(:,3)-params_xy2(5))
%         title(['xy shift, ', num2str(resnorm_xy2)])
% 
%         drawnow                   
%         saveas(gcf,[path2save,int2str(datarunID), '/tiff/ID_',int2str(cellID),'_cones_', int2str(cone1),'_', int2str(cone2),'.tiff'])
%         
%         close(gcf)
%     end
%     
% end
% 
% end


