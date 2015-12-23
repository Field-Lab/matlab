function res=raster(datarun,cell_specification,trigger,varargin)
%
% greschner
% Reverted to Martin's original.  My changes now living in rasterphli.m


% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('start', 0);
    p.addParamValue('stop', []);
    p.addParamValue('stop_trigger_interval', 1);
    p.addParamValue('color', [0 0 0]);
    p.addParamValue('scale', 1);
    p.addParamValue('MarkerSize', 10);
    p.addParamValue('hist', 0);
    p.addParamValue('hist_bin', .01);
    p.addParamValue('hist_color', [1 0 0]);
    p.addParamValue('YLim', []);
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;
    

index=get_cell_indices(datarun,cell_specification);

if isempty(params.stop)
    params.stop=mean(diff(trigger))*params.stop_trigger_interval;
end

res=[];


for i=1:length(trigger),
    h=datarun.spikes{index}-trigger(i);  
    hh=find(h>=params.start & h<=params.stop);
    res=[res; (h(hh)*1000), ones(size(hh))*(length(trigger)-i)];
end

if ~isempty(res)
    if params.hist
        bin=[params.start:params.hist_bin:params.stop];
        t=zeros(length(trigger),length(bin));

        for i=1:length(trigger)
            h=datarun.spikes{index}-trigger(i);  
            t(i,:)=histc(h,bin);
        end
        tt=mean(t,1)/params.hist_bin;

        [AX,H1,H2] = plotyy(res(:,1),res(:,2)*params.scale,(bin+params.hist_bin/2)*1000, tt,'plot');
        set(H1,'LineStyle','.','color',params.color, 'MarkerSize',params.MarkerSize)
        set(H2,'LineStyle','-','color',params.hist_color)
        set(AX,'XLim',[params.start*1000 params.stop*1000],'YColor', [0 0 0]);
        if ~isempty(params.YLim)
            set(AX(2),'YLim',params.YLim);
        end
    else
        plot(res(:,1),res(:,2)*params.scale,'.','color',params.color, 'MarkerSize',params.MarkerSize);
        axis([params.start*1000 params.stop*1000 0 length(trigger)*params.scale]);
    end

end







