
function [mov_reshape,xlim,ylim] = movie_reshape(mov_show,s1,s2,mask)
% s1=40; s2 = 40;
mm = reshape(1:s1*s2,[s1,s2]);
iidx = mm(mask);
[r,c] = find(repelem(reshape(mask,[s1,s2])',1,1));
rmin=min(r);rmax=max(r);cmin=min(c);cmax=max(c);

u=zeros(s1*s2,1);

%mov_reshape=zeros(rmax-rmin+1,cmax-cmin+1,size(mov_show,2));
mov_reshape = zeros(s1,s2,size(mov_show,2));

for itime=1:size(mov_show,2)
u(iidx) = mov_show(:,itime);
xx = reshape(u,[s1,s2]);
%mov_reshape(:,:,itime) = xx(rmin:rmax,cmin:cmax);
mov_reshape(:,:,itime) = xx;
end

xlim=[rmin,rmax];
ylim=[cmin,cmax];
end