function [f g] = regularize_nonsep(p,basepars,f1,g1)

f = f1;
g = g1;

npars_perneuron = length(p)/basepars.Nneurons;

for j=1:basepars.Nneurons

    offset = (j-1)*npars_perneuron + 1;
    stim_idx = offset+1:offset+basepars.n*basepars.Mk;
    
    Kj = p(stim_idx);
    f = f + basepars.lambda_reg/2 * norm(Kj,'fro')^2;
    
    g(stim_idx) = g(stim_idx) + basepars.lambda_reg .* p(stim_idx);
    
end

