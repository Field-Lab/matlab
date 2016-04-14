function [metric,metric_su] = smoothness_metric2(u_spatial_log,mask2)

u_sp_norm = u_spatial_log./repelem(sqrt(sum(u_spatial_log.^2,1)),size(u_spatial_log,1),1);

%% make the shape of filters 2D

sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));

metric_su=[];
%figure;
for isu = 1:size(u_sp_norm,2)
%    subplot(1,size(u_sp_norm,2),isu);
     u_spatial = reshape_vector(u_sp_norm(:,isu),masked_frame,indexedframe);
     u_spatial = u_spatial/max(abs(u_spatial(:)));
%     imagesc(u_spatial);axis image;colormap gray;hold on;
     [r,c] = find(abs(u_spatial) == max(abs(u_spatial(:))));
     u_spatial = u_spatial *sign(u_spatial(r,c));
     
     x=repelem(1:sta_dim2,sta_dim1,1);
     y=repelem([1:sta_dim1]',1,sta_dim2);
     [fitresult, zfit, fiterr, zerr, resnorm, rr] = fmgaussfit(x,y,u_spatial)
     metric_su=[metric_su;resnorm];
     xy=cell(2);xy{1} = x(:);xy{2} = y(:);zz= gaussian2D(fitresult,xy);zz = reshape(zz,[sta_dim1,sta_dim2]);
%     contour(zz,'r');
end

metric = rms(metric_su);
end

function z = gaussian2D(par,xy)
% compute 2D gaussian
z = par(7) + ...
    par(1)*exp(-(((xy{1}-par(5)).*cosd(par(2))+(xy{2}-par(6)).*sind(par(2)))./par(3)).^2-...
    ((-(xy{1}-par(5)).*sind(par(2))+(xy{2}-par(6)).*cosd(par(2)))./par(4)).^2);
end
