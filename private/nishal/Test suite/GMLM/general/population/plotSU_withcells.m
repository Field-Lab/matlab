function [B_use,h,B_use_lr]=plotSU_withcells(K,mask,total_mask_log,Bwt)
K2=gather(K);
s1=size(mask,1); %40,64
s2 = size(mask,2);  %40,32
mm = reshape(1:s1*s2,[s1,s2]);

iidx = mm(mask);
[r,c] = find(repelem(reshape(mask,[s1,s2])',20,20));

ncell = size(total_mask_log,2);
B_use = cell(ncell,1);
for icell=1:ncell
 B = bwboundaries(repelem(reshape(total_mask_log(:,icell),[s1,s2]),20,20));
B_use{icell} = B{1};
end

B_use_lr = cell(ncell,1);
for icell=1:ncell
 B = bwboundaries(repelem(reshape(total_mask_log(:,icell),[s1,s2]),1,1));
B_use_lr{icell} = B{1};
end


cols = distinguishable_colors(ncell+10); cols = cols(5:end,:);

h=[];
%h=figure('Color','w');
nSU = size(K2,2);
for isu=1:size(K2,2)
    subplot(ceil(sqrt(nSU)),ceil(sqrt(nSU)),isu);
    u=zeros(s1*s2,1);
    u(iidx) = K2(:,isu);
    xx = reshape(u,[s1,s2])';
    xx = repelem(xx,20,20);
    
    xx=repelem(xx,1,1,3);

    
    
    imagesc(xx(min(r):max(r),min(c):max(c)));
    % imagesc(xx);
    
    hold on;
    for icell=1:ncell
        ii = convhull(B_use{icell}(:,1),B_use{icell}(:,2),'simplify',true);
    plot(B_use{icell}(ii,1)-min(c)+1,B_use{icell}(ii,2)-min(r)+1,'Color',cols(icell,:),'LineWidth',2);
    end
    
    axis image
    set(gca,'xTick',[]);
    set(gca,'yTick',[]);
    colormap gray
    titleStr ='';
    for icell = 1:ncell
    titleStr = [titleStr,sprintf('%0.02f,',Bwt(isu,icell))];
    end
    title(titleStr);
end

end