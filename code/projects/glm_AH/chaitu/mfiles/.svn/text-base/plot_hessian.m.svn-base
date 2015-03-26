% Plot the gradient values

function out = plot_hessian(x,basepars,stimpars,trainpars)


out = 0;
if (strcmp(basepars.filtermode,'fixfilt'))
    [lg_p cifs] = train_ll3_fixfilt(x,basepars,stimpars,trainpars);
else
    [lg_p cifs] = train_ll3(x,basepars,stimpars,trainpars);
end
    
b_hessval = -trainpars.dt .* sum(cifs);
ps_hess_vals = zeros(basepars.nofilters_postspike,1);
cp_hess_vals = zeros(basepars.nofilters_coupling,1);

% Calculate stimulus filter hess vals

if (strcmp(basepars.filtermode,'sep_raw') || strcmp(basepars.filtermode,'sep_basis'))
    
    [nspace ntime] = get_nspacetime(basepars);

    spcefilt_hess_vals1 = zeros(nspace,1);
    spcefilt_hess_vals2 = zeros(nspace,1);
    tmpfilt_hess_vals1 = zeros(ntime,1);    
    tmpfilt_hess_vals2 = zeros(ntime,1);    
    
     L = linearfilt_grad(x,1,basepars,stimpars,trainpars,1);
     kernel = -trainpars.dt.*block_spikes(cifs',basepars.fac)';
     
     for j=1:nspace
         spcefilt_hess_vals1(j) = dot(L(j,:).^2,kernel);
         spcefilt_hess_vals2(j) = dot(L(nspace+ntime+j,:).^2,kernel);
     end
     
     
     for j=1:ntime
         tmpfilt_hess_vals1(j) = dot(L(nspace+j,:).^2,kernel);
         tmpfilt_hess_vals2(j) = dot(L(nspace+ntime+nspace+j,:).^2,kernel);
     end  
     
     if 0
     H = mult_col(L,kernel)*L';
     offset = get_pars_idx(basepars,1,size(trainpars.D,2),'k');
     offset = offset(1)-1;
     [spce1_idx tmp1_idx spce2_idx tmp2_idx] = get_sep_filt_idces(offset,basepars);
          
     spce1_cond = cond(H(spce1_idx,spce1_idx));
     tmp1_cond = cond(H(tmp1_idx,tmp1_idx));
     spce2_cond = cond(H(spce2_idx,spce2_idx));
     tmp2_cond = cond(H(tmp2_idx,tmp2_idx));
     
     [spce1_cond tmp1_cond spce2_cond tmp2_cond]
     
     fprintf('Condition number of H is %f\n',cond(H));
     
     %H = inv(H);
     %spce1_cond = cond(H(spce1_idx,spce1_idx));
     %tmp1_cond = cond(H(tmp1_idx,tmp1_idx));
     %spce2_cond = cond(H(spce2_idx,spce2_idx));
     %tmp2_cond = cond(H(tmp2_idx,tmp2_idx));
     
     %[spce1_cond tmp1_cond spce2_cond tmp2_cond]     
     
     %fprintf('Condition number of Hinv is %f\n',cond(H));
     end
         
end





% Calculate PS hess vals
for j=1:basepars.nofilters_postspike
    ps_hess_vals(j) = -trainpars.dt*dot(trainpars.psbasisGrad{trainpars.baseneuron_idx(1)}(j,:).^2,cifs);
end


% Calculate CP hess vals
for j=1:basepars.nofilters_coupling
    
    for k=[1:(trainpars.baseneuron_idx-1) trainpars.baseneuron_idx+1:size(trainpars.D,2)]
        cp_hess_vals(j) = cp_hess_vals(j) +abs(trainpars.dt*dot(trainpars.cpbasisGrad{k}(j,:).^2,cifs));
    end
    cp_hess_vals(j) = cp_hess_vals(j)/(size(trainpars.D,2)-1);
end

cla;
bar(1,log(abs(b_hessval)),'b');
offset = 1;

if (strcmp(basepars.filtermode,'sep_basis'))
        for j=1:nspace
            hold on, bar(offset+1,log(abs(spcefilt_hess_vals1(j))),'r');
            offset = offset + 1;
        end
        for j=1:ntime
            hold on, bar(offset+1,log(abs(tmpfilt_hess_vals1(j))),'y');
            offset = offset + 1;
        end        
        for j=1:nspace
            hold on, bar(offset+1,log(abs(spcefilt_hess_vals2(j))),'r');
            offset = offset + 1;
        end        
        for j=1:ntime
            hold on, bar(offset+1,log(abs(tmpfilt_hess_vals2(j))),'y');
            offset = offset + 1;
        end        
    
    
end

if (strcmp(basepars.filtermode,'sep_raw'))
    hold on, bar(offset+1,log(mean(abs(spcefilt_hess_vals1))),'r');
    offset = offset + 1;%basepars.Mk;
    hold on, bar(offset+1,log(mean(abs(tmpfilt_hess_vals1))),'y');
    offset = offset + 1;%basepars.Mk;
    hold on, bar(offset+1,log(mean(abs(spcefilt_hess_vals2))),'r');
    offset = offset + 1;%basepars.Mk;
    hold on, bar(offset+1,log(mean(abs(tmpfilt_hess_vals2))),'y');
    offset = offset + 1;%basepars.Mk;
end
if (basepars.nofilters_postspike > 0)
    hold on, bar(offset+1:offset+basepars.nofilters_postspike,log(abs(ps_hess_vals)),'g');
    offset = offset + basepars.nofilters_postspike;
end
if (size(trainpars.D,2) > 1 && basepars.nofilters_coupling > 0)
    hold on, bar(offset+1:offset+basepars.nofilters_coupling,log(abs(cp_hess_vals)),'c');
end
