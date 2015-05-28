function K = makeToeplitz(templates,Tfind,T)
K=[];


nNeurons = length(templates);
E        = size(templates{1},1);


for n = 1:nNeurons

    for t = 1:length(Tfind)
        
        ActionPotential            = makeActionPotential(n,Tfind(t),templates,T);
        K(:,t+(n-1)*length(Tfind)) = reshape(ActionPotential',E*T,1);
end

end
    