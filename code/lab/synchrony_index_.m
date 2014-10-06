function sync=synchrony_index_(sp1,sp2, varargin) 
% greschner


% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('bin', 10);% ms
    p.addParamValue('duration', []);% ms
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;
 
    
if isempty(params.duration)
    params.duration=max([sp1;sp2]);
end

bin=0:params.bin/1000:params.duration-params.bin/1000;

spiketrain=zeros(length(bin),2);
spiketrain(:,1)=histc(sp1,bin);
spiketrain(:,2)=histc(sp2,bin);

temp=corrcoef(spiketrain);
sync=temp(1,2);
         
















    
    