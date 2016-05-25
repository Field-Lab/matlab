% It returns function, gradient (and Hessian?) value of the log likelihood
% of the data with respect to the base firing rate, MAGNITUDE of the stimulus filter with fixed shape k, and postspike/coupling filters coeffs.
% This function can be used by fmincon

function [f g Hinfo] = ll_func2_fixfilt_hessmult(p,k,kx_base,basepars,stimpars,trainpars)
% Calculate the log likelihood of the data


nfullpars_perneuron = 1+basepars.n*basepars.Mk+basepars.nofilters_postspike + (basepars.Nneurons-1)*basepars.nofilters_coupling;
npartpars_perneuron = 2+basepars.nofilters_postspike + (basepars.Nneurons-1)*basepars.nofilters_coupling;
pfull = magtofull(p,basepars,k);
        
[logprob cifs] = train_ll3(pfull,basepars,stimpars,trainpars,p(2).*kx_base);
f = full(-logprob); % remember to negate - we want to maximize log likelihood
% Calculate gradient
g = zeros(size(p));
kernel1 = full((1.*trainpars.D) - trainpars.dt.*cifs);

if (nargout > 2)
   
    Hinfo.kernel2 = full(-trainpars.dt.*cifs);
    %k2blocked = block_spikes(kernel2',basepars.fac)';
    
    %Hinfo.kernel1 = kernel1;
    Hinfo.kx = kx_base;

end

1;
for j=1:basepars.Nneurons
    offset = npartpars_perneuron*(j-1);
    g(offset+1) = sum(kernel1(:,j));
    g(offset+2) = dot(kx_base(:,j),block_spikes(kernel1(:,j)',basepars.fac));
    g(offset+3:offset+2+basepars.nofilters_postspike) = trainpars.psbasisGrad(:,:,j)*kernel1(:,j);
    
    if (basepars.Nneurons > 1)
        cpGrad = zeros(basepars.nofilters_coupling*(basepars.Nneurons-1),size(kernel1,1));
        for n2 = [1:j-1 j+1:basepars.Nneurons] % neuron n is the one being AFFECTED
            n2rev = n2;
            if (n2 > j)
                n2rev = n2rev-1;
            end
            cpGrad((n2rev-1)*basepars.nofilters_coupling+1:n2rev*basepars.nofilters_coupling,:) = trainpars.cpbasisGrad(:,:,n2);
        end
        
        g(offset+2+basepars.nofilters_postspike+1:j*npartpars_perneruon) = cpGrad*kernel1;
        
    end
    
    

    
    
    
end

g = -g;

fprintf('f=%f    gnorm=%f\n',f,norm(g));





