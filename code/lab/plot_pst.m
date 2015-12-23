function [psth]=plot_pst(datarun, cell_id, tr, varargin) 
%greschner

% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('psth', 1);
    p.addParamValue('dt', .001);
    p.addParamValue('raster_color', 'k');
    p.addParamValue('marker_size', 6);
    p.addParamValue('hist_color', 'r');
    p.addParamValue('start', 0);
    p.addParamValue('stop', []);
    p.addParamValue('yscale', 1);
    p.addParamValue('gauss_smooth', 0);
    p.parse(varargin{:});
    params = p.Results;

    
if isempty(params.stop) 
    params.stop=mean(diff(tr));
end
 
index=get_cell_indices(datarun,cell_id);
spikes=datarun.spikes{index};
  
psth_r=[];
for i=1:length(tr),
    h=spikes-tr(i);  
    hh=find(h>=params.start & h<=params.stop);
    psth_r=[psth_r; (h(hh)*1000)', repmat(length(tr)-i,[length(hh),1]);];
end

if ~isempty(psth_r)
    plot(psth_r(:,1),psth_r(:,2)*y_scale,[raster_color '.'],'MarkerSize',params.marker_size);
    axis([mmin*1000 mmax*1000 0 length(tr)*y_scale]);
end


if params.psth
    
    bin=[params.start:params.dt:params.stop]
    t=zeros(length(tr),length(bin));
    for i=1:length(tr)
       t(i,:)=histc(spikes,tr(i)+bin);  
    end
    psth=sum(t,1);
    
    if params.gauss_smooth~=0
        psth=gauss_smooth(psth,params.gauss_smooth);
    end
    
    tt=psth/max(psth)*y_scale*length(tr); 
    hold on
    plot(bin,tt,params.hist_color);

else
    psth=[];
end








   
    






























    
    