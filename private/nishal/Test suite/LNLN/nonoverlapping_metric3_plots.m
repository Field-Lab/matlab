function  [fits,clog]= nonoverlapping_metric3_plots(u_spatial_log,mask2,sta)

u_sp_norm = u_spatial_log./repelem(sqrt(sum(u_spatial_log.^2,1)),size(u_spatial_log,1),1);

%% make the shape of filters 2D

sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));

metric_su=[];
%figure;
nSU = size(u_sp_norm,2);
lw=1.5;
[r,c] = find(mask2>0);
u_sp_log = cell(nSU,1);
    x_coord = [min(r)-2:max(r)+2];
    y_coord = [(min(c)-2):(max(c)+2)];
    
xxsta =-repelem(sta.*mask2,20,20);
B = bwboundaries(abs(xxsta)>0.2*max(abs(xxsta(:))));
hullidx = convhull(B{1}(:,1),B{1}(:,2),'simplify',true);
xxsta = xxsta/max(abs(xxsta(:)));

figure;
subplot(1,nSU+2,1);
imagesc((1-repelem(-sta/max(abs(sta(:))),1,1,3))/2);axis image;hold on;
ylim([min(r)-2,max(r)+2]);
xlim([min(c)-2,max(c)+2]);hold on;
plot((B{1}(hullidx,2)/20)+0.5,(B{1}(hullidx,1)/20)+0.5,'LineWidth',lw);set(gca,'xTick',[]);set(gca,'yTick',[]);

for isu = 1:nSU

     u_spatial = reshape_vector(u_sp_norm(:,isu),masked_frame,indexedframe);
     u_spatial = u_spatial/max(abs(u_spatial(:)));
     u_sp_log{isu} = u_spatial;
     
     
     x=repelem(1:sta_dim2,sta_dim1,1);
     y=repelem([1:sta_dim1]',1,sta_dim2);
     [fitresult, zfit, fiterr, zerr, resnorm, rr] = fmgaussfit(x,y,u_spatial)
     metric_su=[metric_su;resnorm];
     xy=cell(2);xy{1} = x(:);xy{2} = y(:);zz= gaussian2D(fitresult,xy); zz = reshape(zz,[sta_dim1,sta_dim2]);
     
     fits(isu).fitresult = fitresult;  
         
     subplot(1,size(u_sp_norm,2)+2,isu+1);
     imagesc((1 - repelem(u_spatial,1,1,3))/2);axis image;colormap gray;hold on;
     plot((B{1}(hullidx,2)/20)+0.5,(B{1}(hullidx,1)/20)+0.5,'LineWidth',lw);hold on;
     contour(zz,'r');hold on
     ylim([min(r)-2,max(r)+2]);
     xlim([min(c)-2,max(c)+2]);
     caxis([0,1]);set(gca,'xTick',[]);set(gca,'yTick',[]);
     
     
end

subplot(1,nSU+2,nSU+2);
xxx = [0:0.1:sta_dim2];
yyy = [0:0.1:sta_dim1];
imge=zeros(numel(yyy),numel(xxx));
x=repelem(0:0.1:sta_dim2,numel(yyy),1);
y=repelem([0:0.1:sta_dim1]',1,numel(xxx));

xy=cell(2,1);xy{1} = x(:);xy{2} = y(:);
cols = distinguishable_colors(nSU);
for isu=1:nSU
    zz= gaussian2D(fits(isu).fitresult ,xy);
    iidx = zz>normpdf(1,0,1)*max(zz);
    % plot(xy{1}(iidx),xy{2}(iidx),'.','Color',cols(isu,:));
    %    hold on;
    imge=imge+ (reshape(iidx,[numel(yyy),numel(xxx)]));
    zz=zz/sum(zz);
    fits(isu).mean2 = [zz'*xy{1}, zz'*xy{2}]
    fits(isu).mean = [mean(xy{1}(iidx)), mean(xy{2}(iidx))];
end
imagesc(xxx,yyy,(1-repelem(imge/max(abs(imge(:))),1,1,3))/2);axis image;hold on;caxis([0,1]);
hold on;
plot((B{1}(hullidx,2)/20)+0.5,(B{1}(hullidx,1)/20)+0.5,'LineWidth',lw);
set(gca,'xTick',[]);set(gca,'yTick',[]);

ylim([min(r)-2,max(r)+2]);
xlim([min(c)-2,max(c)+2]);
     
for isu=1:nSU
    plot(fits(isu).mean(1),fits(isu).mean(2),'r.','MarkerSize',20);hold on;
end


 ylim([(min(r)-2),(max(r)+2)]);
 xlim();
     
 %% relative distance between centers, normalize it 
 center_log=[];othercenter_log=[];
 for isu=1:nSU
     for jsu=1:nSU
         if(isu == jsu);continue;end
         xy=cell(2,1);
         xy{1} = fits(isu).mean(1);
         xy{2} = fits(isu).mean(2);
         center = gaussian2D(fits(isu).fitresult,xy);
         
         xy=cell(2,1);
         xy{1} = fits(jsu).mean(1);
         xy{2} = fits(jsu).mean(2);
         othercenter = gaussian2D(fits(isu).fitresult,xy);
         
         center_log= [center_log;center];
         othercenter_log= [othercenter_log;othercenter];
     end
 end
clog.cc_ratio =othercenter_log./center_log;
end


function z = gaussian2D(par,xy)
% compute 2D gaussian
z = par(7) + ...
    par(1)*exp(-(((xy{1}-par(5)).*cosd(par(2))+(xy{2}-par(6)).*sind(par(2)))./par(3)).^2-...
    ((-(xy{1}-par(5)).*sind(par(2))+(xy{2}-par(6)).*cosd(par(2)))./par(4)).^2);
end


