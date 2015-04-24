function [sta nsp waveforms] = fast_sta(X,sp_times,m,dt,rotateFlag)

if (~exist('rotateFlag','var'))
    rotateFlag = 0;
end

[n T] = size(X);
sta = zeros(n,m);

tvec = 0:dt:(T-1)*dt;
nsp = 0;
1;

if (nargout > 2)
    waveforms = zeros(m,length(sp_times));
end

for i=1:length(sp_times)
    
    idx = find(tvec > sp_times(i),1);
    if (isempty(idx) || idx < m)
        continue;
    end
    
    sta = sta + X(:,idx:-1:idx-m+1);
    
    if (nargout > 2)
        waveforms(:,i) = X(:,idx:-1:idx-m+1);
    end
    
    nsp = nsp+1;
end
fprintf('Fast STA: used %d spikes.\n',nsp);
sta = sta./nsp;

MAX_SIZE = 10^2;
1;
if (rotateFlag) % rotate by the inverse autocorrelation matrix - this is potentially VERY expensive
    Xdes = stim_design(X,m);
    sta = reshape((Xdes(:,1:MAX_SIZE)*Xdes(:,1:MAX_SIZE)')\sta(:),n,m);
end