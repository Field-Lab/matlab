idx = nchoosek(1:k,m-n);

for i = 1:length(idx)
    AD = [A D(:,idx(i,:))];
    if rank(AD)==m
        break;
    end
end

% B is the equalizer
% Drej is the set of rejected disturbance patterns
BE = inv(AD);
B = BE(1:n,:);
Drej = AD(:,(n+1):m);

% Checks:
% size(AD) = mxn;
% rank(AD) = m;
% B*A = I;
% B*Drej = 0;