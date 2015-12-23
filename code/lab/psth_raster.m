function psth_r=psth_raster(mmin,mmax,n,tr,Color,y_scale)
% plots pst spike raster
% psth_raster(mmin,mmax,spikes,trigger,'k')

psth_r=[];

%if nargin==5
%    for a=1:nint-1;
%        plot([mmin*1000 mmax*1000],[a*length(tr)/nint a*length(tr)/nint],'Color',[.9 .9 .9]);
%        hold on;
%    end  
%end

if nargin==4
    Color='k';
end
Color=[Color '.'];

if nargin<6
    y_scale=1;
else
    y_scale=y_scale/length(tr);
end

for i=1:length(tr),
        h=n-tr(i);  
        hh=find(h>=mmin & h<=mmax);
        psth_r=[psth_r; (h(hh)*1000)', repmat(length(tr)-i,[length(hh),1]);];
end

if ~isempty(psth_r)
    plot(psth_r(:,1),psth_r(:,2)*y_scale,Color);%, 'MarkerSize',10
    axis([mmin*1000 mmax*1000 0 length(tr)*y_scale]);
end
hold off;