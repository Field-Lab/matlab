function [ccf time ccf_raw ccf_shuffle]=correlation_shift_shuffle(sp1, sp2, triggers, start ,stop, varargin) 
% shifted corrcoeff eg crosscorrelation function
%
%
%
% greschner

% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('dt', .004);
    p.addParamValue('offset', .15);
    p.addParamValue('shift', []);
    p.parse(varargin{:});
    params = p.Results;

    
if isempty(params.shift)
    params.shift=params.dt;
end    
    

n=round(params.offset/params.shift);
n=-n:n;
time=n*params.shift;
ccf_raw=zeros(size(time));
ccf_shuffle=zeros(size(time));

spiketrain=zeros((stop-start)/params.dt,2*length(triggers));
tt1= logical(eye(length(triggers)));

for j=1:length(triggers)
    spiketrain(:,j)=histc(sp1,triggers(j)+start:params.dt:triggers(j)+stop-params.dt); 
end

for jj=1:length(time)
    for j=1:length(triggers)
        spiketrain(:,length(triggers)+j)=histc(sp2-time(jj),triggers(j)+start:params.dt:triggers(j)+stop-params.dt); 
    end
    tem=corrcoef(spiketrain);
    temp=tem(1:length(triggers),length(triggers)+1:end);

    ccf_raw(jj)=mean(temp(tt1)); 
    ccf_shuffle(jj)=mean(temp(~tt1)); 
end

ccf=ccf_raw-ccf_shuffle;  






    

















    
    