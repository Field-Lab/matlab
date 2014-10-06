function [ccf, time]=get_xcorr(data, cell_id1, cell_id2, varargin)
% greschner

% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('dt', .002);
    p.addParamValue('offset', .75);
    p.addParamValue('stop', []);
    p.addParamValue('start', []);
    p.parse(varargin{:});
    params = p.Results;

    
index1=get_cell_indices(data,cell_id1);
index2=get_cell_indices(data,cell_id2); 
if isempty(index1) | isempty(index1),
    error('get_xcorr: cell_id does not exit');
end

sp1=data.spikes{index1};
sp2=data.spikes{index2};
         
if isempty(params.start)
    params.start=min([sp1(1) sp2(1)]);
end    
if isempty(params.stop)
    params.stop=max([sp1(end) sp2(end)]);
end
    

hsp1=histc(sp1,[params.start:params.dt:params.stop]);    
hsp2=histc(sp2,[params.start:params.dt:params.stop]);

[ccf, lags]=xcorr(hsp2,hsp1,round(params.offset/params.dt),'coeff');
time=lags*params.dt;

