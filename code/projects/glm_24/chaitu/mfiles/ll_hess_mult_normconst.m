function W = ll_hess_mult_normconst(x,lambda,Y,basepars,stimpars,trainpars)

% Function that indirectly compute H*Y where H is the hessian of log
% likelihood

% Hinfo.kernel2 is the kernel (Tx1) vector (weighting vector)


% The trick is that H = A*K*A' where K is diagonal matrix
% Computing H*Y conventionally takes:
%    O(nT) for computing Hinfo*K;
%    O(n^2*T) for computing (Hinfo*K) * Hinfo'
%    O(n^2*b) for computing H*Y
% so the dominant term is O(n^2*T)

% Indirectly it takes:
%   O(n*T) for computing Hinfo*K
%   O(n*T*b) for computing Hinfo'*Y
%   O(n*T*b) for computing (Hinfo*K) * (Hinfo'*Y)
% so the dominant term is O(n*T*b) which is much less if n is large


% It will be assumed that ther diagonal entries of the K matrix will be
% stored in Hinfo.kernel2

%fprintf('I am being called\n');
[lgp lcifs cifs] = train_ll3(x,basepars,stimpars,trainpars);
fprime = basepars.Nprime(lcifs);
fdprime = basepars.Ndoubleprime(lcifs);
Hinfo.kernel2 = (fdprime.*cifs - fprime.^2)./(cifs.^2).*full(double(trainpars.D(:,trainpars.baseneuron_idx))) - fdprime.*trainpars.dt;
W = ll_hess_mult(Hinfo,Y,basepars,stimpars,trainpars); % use the regular hess mult fn

nneurons = basepars.Nneurons;
Neff = size(trainpars.D,2);
sqflag = (isfield(basepars,'XsqK') && basepars.XsqK);

% Add the lagrangian terms

if (~isfield(basepars,'normconstr_p'))
    basepars.normconstr_p = 2;
end
p = basepars.normconstr_p;

switch(basepars.filtermode)
    
    case 'nonsep'
        
        for j=1:nneurons
            stim_idx  = get_pars_idx(basepars,j,Neff,'k');
            W(stim_idx,:) = W(stim_idx,:) + lambda.ineqnonlin(j).*((p-1)*abs(x(stim_idx)).^(p-2)).*Y(stim_idx,:);
            
            if (sqflag)
                sqstim_idx = get_pars_idx(basepars,j,Neff,'ksq');
                W(sqstim_idx,:) = W(sqstim_idx,:) + lambda.ineqnonlin(nneurons+j).*((p-1)*abs(x(sqstim_idx)).^(p-2)).*Y(sqstim_idx,:);
            end
            
        end
        
    case {'sep_raw','sep_basis'}
        for j=1:nneurons
            stim_idx = get_pars_idx(basepars,j,Neff,'k');
            [stim_idx1 blah stim_idx2 blah] = get_sep_filt_idces(stim_idx(1)-1,basepars);

            W(stim_idx1,:) = W(stim_idx1,:) + lambda.ineqnonlin(2*j-1).*((p-1)*abs(x(stim_idx1)).^(p-2)).*Y(stim_idx1,:);
            W(stim_idx2,:) = W(stim_idx2,:) + lambda.ineqnonlin(2*j).*((p-1)*abs(x(stim_idx2)).^(p-2)).*Y(stim_idx2,:);
            
            if (sqflag)
               stim_idx = get_pars_idx(basepars,j,Neff,'ksq');
               [stim_idx1 blah stim_idx2 blah] = get_sep_filt_idces(stim_idx(1)-1,basepars);

                W(stim_idx1,:) = W(stim_idx1,:) + lambda.ineqnonlin(2*nneurons+2*j-1).*((p-1)*abs(x(stim_idx1)).^(p-2)).*Y(stim_idx1,:);
                W(stim_idx2,:) = W(stim_idx2,:) + lambda.ineqnonlin(2*nneurons+2*j).*((p-1)*abs(x(stim_idx2)).^(p-2)).*Y(stim_idx2,:);

            end
            
        end

        
        
end
    
    
    
    
    

