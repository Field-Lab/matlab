function [f g H] = regularize_sep_whole(p,basepars,f1,g1,H1)

if (~strcmp(basepars.filtermode,'sep_raw'))
    fprintf('ERROR: cannot regularize  separable filters in nonsep mode!\n');
    return;
end

f =f1;
g = g1;

if (nargout >= 3)
    H = H1;
end

npars_perneuron = length(p)/basepars.Nneurons;

for j=1:basepars.Nneurons
    
    offset = (j-1)*npars_perneuron+1;

    % Stimulus filter component indices
    [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(offset,basepars);

    % Regularize the whole filter in nonsep mode (will do more sophisticated options later)
    Kj = p(s1_idx)*p(t1_idx)' + p(s2_idx)*p(t2_idx)';
    f = f + basepars.lambda_reg/2 * norm(Kj,'fro')^2;
    g(s1_idx) = g(s1_idx) + (basepars.lambda_reg*norm(p(t1_idx))^2).*p(s1_idx);
    g(t1_idx) = g(t1_idx) + (basepars.lambda_reg*norm(p(s1_idx))^2).*p(t1_idx);
    g(s2_idx) = g(s2_idx) + (basepars.lambda_reg*norm(p(t2_idx))^2).*p(s2_idx);
    g(t2_idx) = g(t2_idx) + (basepars.lambda_reg*norm(p(s2_idx))^2).*p(t2_idx);
    
    if (nargout >= 3)
        
        % Diagonal terms
        H(s1_idx,s1_idx) = H(s1_idx,s1_idx) + (basepars.lambda_reg*norm(p(t1_idx))^2).*eye(length(s1_idx));
        H(t1_idx,t1_idx) = H(t1_idx,t1_idx) + (basepars.lambda_reg*norm(p(s1_idx))^2).*eye(length(t1_idx));
        H(s2_idx,s2_idx) = H(s2_idx,s2_idx) + (basepars.lambda_reg*norm(p(t2_idx))^2).*eye(length(s2_idx));
        H(t2_idx,t2_idx) = H(t2_idx,t2_idx) + (basepars.lambda_reg*norm(p(s2_idx))^2).*eye(length(t2_idx));
        
        % Cross-terms
        H(s1_idx,t1_idx) = H(s1_idx,t1_idx) + (2*basepars.lambda_reg).*(p(s1_idx)*p(t1_idx)');
        H(t1_idx,s1_idx) = H(t1_idx,s1_idx) + (2*basepars.lambda_reg).*(p(t1_idx)*p(s1_idx)');
        H(s2_idx,t2_idx) = H(s2_idx,t2_idx) + (2*basepars.lambda_reg).*(p(s2_idx)*p(t2_idx)');
        H(t2_idx,s2_idx) = H(t2_idx,s2_idx) + (2*basepars.lambda_reg).*(p(t2_idx)*p(s2_idx)');
        
    end
    
end