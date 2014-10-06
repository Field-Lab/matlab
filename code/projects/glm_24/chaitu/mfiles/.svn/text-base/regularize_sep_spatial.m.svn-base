function [f g H] = regularize_sep_spatial(p,basepars,f1,g1,H1)

if (~strcmp(basepars.filtermode,'sep_raw'))
    fprintf('ERROR: cannot regularize  spatial filters in nonsep mode!\n');
    return;
end

npars_perneuron = length(p)/basepars.Nneurons;

f = f1;
g = g1;

if (nargout >= 3)
    H = H1;
end

for j=1:basepars.Nneurons
    offset = (j-1)*npars_perneuron+1;
    % Stimulus filter component indices
    [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(offset,basepars);

    % Regularize only the spatial filters!
    f = f + basepars.lambda_reg/2 * (norm(p(s1_idx),'fro')^2 + norm(p(s2_idx),'fro')^2);
    g(s1_idx) = g(s1_idx) + basepars.lambda_reg.*p(s1_idx);
    g(s2_idx) = g(s2_idx) + basepars.lambda_reg.*p(s2_idx);

    if (nargout >= 3)
       % Diagonal terms
       H(s1_idx,s1_idx) = H(s1_idx,s1_idx) + basepars.lambda_reg.*eye(length(s1_idx));
       H(s2_idx,s2_idx) = H(s2_idx,s2_idx) + basepars.lambda_reg.*eye(length(s2_idx));
    end
    
end