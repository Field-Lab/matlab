function plot_mat(m,bin,shift)
if nargin<2
    bin=1;
end
if nargin<3
    shift=[0 0];
end


mt=ones(length(m(:,1))+1,length(m(1,:))+1)*m(1,1);
mt(1:length(m(:,1)),1:length(m(1,:)))=(m);
%pcolor (mt),  shading flat, pbaspect([1 1 1]);
%pcolor (mt),  shading flat;

%h=pcolor (mt);
h=pcolor ([.5*bin+shift(1)*bin:bin:length(mt(1,:))*bin-.5*bin+shift(1)*bin],[.5*bin+shift(2)*bin:bin:length(mt(:,1))*bin-.5*bin+shift(2)*bin],double(mt));
%shading flat
set(h,'LineStyle','none');
%set(gca,'YDir','reverse');
%colormap(gray)
%set(h,'edgecolor','none')