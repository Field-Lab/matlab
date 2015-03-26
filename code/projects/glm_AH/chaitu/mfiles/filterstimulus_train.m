function [kx kxsq] = filterstimulus_train(p,basepars,stimpars,trainpars,fl_init)

%fprintf('\n\n HEY!! filterstimulus_train.m called \n\n')
% FM-DEPENDENT (edoi, 2012-01)

Neff = size(trainpars.D,2);
if (strcmp(basepars.filtermode,'fixrk2') && exist('fl_init','var'))
   kx = cell(2,1);
else
   kx = zeros(size(stimpars.x,2),basepars.Nneurons);
end
if (nargout>1)
    kxsq = zeros(size(kx));
end

T = size(stimpars.x,2);

if (~isfield(basepars,'meta_filt_mode') && strcmp(basepars.filtermode,'stark1'))
   basepars.meta_filt_mode = 'starkn';
end


for j=1:basepars.Nneurons
    
    if (iscell(basepars.crop_idx))
        crop_idx_j = basepars.crop_idx{j};
    else
        crop_idx_j = basepars.crop_idx(:,j);
    end
    
    kidx = get_pars_idx(basepars,j,Neff,'k');
    if (nargout>1)
        ksqidx = get_pars_idx(basepars,j,Neff,'ksq');
    end
    
    
    if ( strcmp(basepars.meta_filt_mode,'indrkn') || strcmp(basepars.meta_filt_mode,'stst') )
       for rnk = 1:basepars.stim_n
          if rnk == 1
             Kfilt = p(2)*basepars.K{1};
          else
             Kfilt = Kfilt + p(rnk+1)*basepars.K{rnk};
          end
       end
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi

    elseif strcmp(basepars.meta_filt_mode,'jit2')
       
       fprintf('### do not need to recompute filtered stimulus\n\n')
       %Ks = p(2:61); Ks = Ks(:);
       %Kfilt = Ks*basepars.Kt;
       %kx = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,60^2,T,basepars.padval),1);
       
    elseif strcmp(basepars.meta_filt_mode,'lcs0')
       Kfilt = p(2)*basepars.K;
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       %-- additional non-linear terms
       tc = basepars.lc0_tc; 
       % 2. convolve with the time course
       K_slen = 15;
       n_sti_fr = size(stimpars.x,2);
       tmp = fastconv(stimpars.x(crop_idx_j,:),repmat(tc,K_slen^2,1),K_slen^2,n_sti_fr,basepars.padval);
       c_reg = zeros(15);
       c_reg(6:10,6:10) = 1;
       c_reg = c_reg(:);
       tmp = tmp(c_reg==1,:);
       % 3. zero-mean
       tmp = tmp - repmat(mean(tmp,2),1,n_sti_fr);
       % 4. 
       tmp_p = zeros(size(tmp));
       tmp_n = zeros(size(tmp));
       tmp_p(tmp>0) = tmp(tmp>0);
       tmp_n(tmp<0) = tmp(tmp<0);
       for k = 1:25
          kx(:,j) = p(2+k)*tmp_p(k,:)'+kx(:,j);
          kx(:,j) = p(2+25+k)*tmp_n(k,:)'+kx(:,j);
       end
       
    elseif strcmp(basepars.meta_filt_mode,'abp')
       Kfilt = p(2)*basepars.K;
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       %-- additional non-linear terms
       %fs_stc = cell(25,1);
       fs_stc = cell(49,1);
       for ck = 1:49 %25
          fs_stc{ck} = stimpars.dt*sum(fastconv(stimpars.x(crop_idx_j,:),basepars.nl_filt{ck},15^2,T,basepars.padval),1);
       end
       ckpn = 0;
       for ck = 1:49 %25
          tmp   = fs_stc{ck};
          tmp_p = zeros(size(tmp));
          tmp_n = zeros(size(tmp));
          tmp_p(tmp>0) = tmp(tmp>0);
          tmp_n(tmp<0) = tmp(tmp<0);
          ckpn = ckpn+1;
          kx(:,j) = p(2+ckpn)*tmp_p.^2'+kx(:,j);
          ckpn = ckpn+1;
          kx(:,j) = p(2+ckpn)*tmp_n.^2'+kx(:,j);
       end
       clear tmp*

    elseif strcmp(basepars.meta_filt_mode,'stc2r')
       Kfilt = p(2)*basepars.K;
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       %-- additional non-linear terms
       n_sti_fr = size(kx,1);
       for k = 1:2
          %tmp = stimpars.dt * sum(fastconv(stimpars.x(crop_idx_j,:),basepars.nl_filt,15^2,n_sti_fr,basepars.padval),1);
          % NOTE: basepars.nl_filt ... bug??
          tmp = stimpars.dt * sum(fastconv(stimpars.x(crop_idx_j,:),basepars.nl_filt{k},15^2,n_sti_fr,basepars.padval),1);
          tmp = tmp - mean(tmp);
          tmp_p = zeros(size(tmp));
          tmp_n = zeros(size(tmp));
          tmp_p(tmp>0) = tmp(tmp>0);
          tmp_n(tmp<0) = tmp(tmp<0);
          kx(:,j) = p(2+1)*tmp_p.^2'+kx(:,j);
          kx(:,j) = p(2+2)*tmp_n.^2'+kx(:,j);
          fprintf('pos/neg separate channe... inconsistent with the ctrl_vglm_fit_ALL_dat4.m, need to debug. \n')
          pause
       end
       clear tmp*
       
    elseif strcmp(basepars.meta_filt_mode,'stc6rs')
       Kfilt = p(2)*basepars.K;
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       %-- additional non-linear terms
       n_sti_fr = size(kx,1);
       c = 0;
       for k = 1:6
          tmp = stimpars.dt * sum(fastconv(stimpars.x(crop_idx_j,:),basepars.nl_filt{k},15^2,n_sti_fr,basepars.padval),1);
          tmp = tmp - mean(tmp);
          tmp_p = zeros(size(tmp));
          tmp_n = zeros(size(tmp));
          tmp_p(tmp>0) = tmp(tmp>0);
          tmp_n(tmp<0) = tmp(tmp<0);
          c = c+1;
          kx(:,j) = p(2+c)*tmp_p.^2'+kx(:,j);
          c = c+1;
          kx(:,j) = p(2+c)*tmp_n.^2'+kx(:,j);
       end
       clear tmp*
       
    elseif strcmp(basepars.meta_filt_mode,'lfst')
       Kfilt = p(2)*basepars.K;
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       %-- additional non-linear terms
       tmp = stimpars.dt * sum(fastconv(stimpars.x(crop_idx_j,:),basepars.nl_filt,15^2,T,basepars.padval),1);
       tmp = tmp - mean(tmp);
       tmp_p = zeros(size(tmp));
       tmp_n = zeros(size(tmp));
       tmp_p(tmp>0) = tmp(tmp>0);
       tmp_n(tmp<0) = tmp(tmp<0);
       tmp_p = basepars.rsc_p*tmp_p.^2;
       tmp_n = basepars.rsc_n*tmp_n.^2;
       % temporal convolution of nonlinear signal ***
       n_tap = basepars.n_nl/2;
       mi = basepars.lfst.tap.mi;
       % pl = basepars.lfst.tap.pl;
       c = 0;
       for k = 1:n_tap
          if k <= mi
             tmp = zeros(1,T);
             tmp(1:end-mi+k-1) = tmp_p(1+mi-k+1:end);
             kx(:,j) = p(2+k)*tmp(:)+kx(:,j);
             tmp = zeros(1,T);
             tmp(1:end-mi+k-1) = tmp_n(1+mi-k+1:end);
             kx(:,j) = p(2+n_tap+k)*tmp(:)+kx(:,j);
          elseif k == mi+1
             kx(:,j) = p(2+k)      *tmp_p+kx(:,j);
             kx(:,j) = p(2+n_tap+k)*tmp_n+kx(:,j);
          else
             c = c+1;
             tmp = zeros(1,T);
             tmp(1+c:end) = tmp_p(1:end-c);
             kx(:,j) = p(2+k)*tmp(:)+kx(:,j);             
             c = c+1;
             tmp = zeros(1,T);
             tmp(1+c:end) = tmp_n(1:end-c);
             kx(:,j) = p(2+n_tap+k)*tmp(:)+kx(:,j);             
          end
       end
       clear tmp

    elseif strcmp(basepars.meta_filt_mode,'lfst0')
       Kfilt = p(2)*basepars.K;
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       %-- additional non-linear terms
       n_sti_fr = size(kx,1);
       tmp = stimpars.dt * sum(fastconv(stimpars.x(crop_idx_j,:),basepars.nl_filt,15^2,n_sti_fr,basepars.padval),1);
       %tmp = tmp*100;
       tmp = tmp - mean(tmp);
       tmp_p = zeros(size(tmp));
       tmp_n = zeros(size(tmp));
       tmp_p(tmp>0) = tmp(tmp>0);
       tmp_n(tmp<0) = tmp(tmp<0);
       tmp_p = basepars.rsc_p*tmp_p.^2;
       tmp_n = basepars.rsc_n*tmp_n.^2;

       n_tap = basepars.n_nl/2;
       mi = basepars.lfst.tap.mi;
       c = 0;
       for k = 1:n_tap
          if k <= mi
             tmp = zeros(1,T);
             tmp(1:end-mi+k-1) = tmp_p(1+mi-k+1:end);
             kx(:,j) = p(2+k)*tmp(:)+kx(:,j);
             tmp = zeros(1,T);
             tmp(1:end-mi+k-1) = tmp_n(1+mi-k+1:end);
             kx(:,j) = p(2+n_tap+k)*tmp(:)+kx(:,j);
          elseif k == mi+1
             kx(:,j) = p(2+k)      *tmp_p(:)+kx(:,j);
             kx(:,j) = p(2+n_tap+k)*tmp_n(:)+kx(:,j);
          else
             c = c+1;
             tmp = zeros(1,T);
             tmp(1+c:end) = tmp_p(1:end-c);
             kx(:,j) = p(2+k)*tmp(:)+kx(:,j);             
             c = c+1;
             tmp = zeros(1,T);
             tmp(1+c:end) = tmp_n(1:end-c);
             kx(:,j) = p(2+n_tap+k)*tmp(:)+kx(:,j);             
          end
       end
       
    elseif strcmp(basepars.meta_filt_mode,'lfs')
       Kfilt = p(2)*basepars.K;
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       %-- additional non-linear terms
       n_sti_fr = size(kx,1);
       tmp = stimpars.dt * sum(fastconv(stimpars.x(crop_idx_j,:),basepars.nl_filt,15^2,n_sti_fr,basepars.padval),1);
       %tmp = tmp*100;
       tmp = tmp - mean(tmp);
       tmp_p = zeros(size(tmp));
       tmp_n = zeros(size(tmp));
       tmp_p(tmp>0) = tmp(tmp>0);
       tmp_n(tmp<0) = tmp(tmp<0);
       tmp_p = tmp_p.^2;
       tmp_n = tmp_n.^2;
       kx(:,j) = basepars.rsc_p*p(2+1)*tmp_p'+kx(:,j);
       kx(:,j) = basepars.rsc_n*p(2+2)*tmp_n'+kx(:,j);
       clear tmp

    elseif strcmp(basepars.meta_filt_mode,'lcs1')
       Kfilt = p(2)*basepars.K;
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       %-- additional non-linear terms
       tc = basepars.lc0_tc; 
       % 2. convolve with the time course
       K_slen = 15;
       n_sti_fr = size(stimpars.x,2);
       tmp = fastconv(stimpars.x(crop_idx_j,:),repmat(tc,K_slen^2,1),K_slen^2,n_sti_fr,basepars.padval);
       c_reg = zeros(15);
       c_reg(6:10,6:10) = 1;
       c_reg = c_reg(:);
       tmp = tmp(c_reg==1,:);
       % 3. zero-mean
       tmp = tmp - repmat(mean(tmp,2),1,n_sti_fr);
       % 4. 
       tmp_p = zeros(size(tmp));
       tmp_n = zeros(size(tmp));
       tmp_p(tmp>0) = tmp(tmp>0);
       tmp_n(tmp<0) = tmp(tmp<0);
       tmp_p = tmp_p.^2;
       tmp_n = tmp_n.^2;
       for k = 1:25
          kx(:,j) = basepars.resc(k)*p(2+k)*tmp_p(k,:)'+kx(:,j);
          kx(:,j) = basepars.resc(k+25)*p(2+25+k)*tmp_n(k,:)'+kx(:,j);
       end
       
    elseif strcmp(basepars.meta_filt_mode,'lcs1t0')
       Kfilt = p(2)*basepars.K;
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       %-- additional non-linear terms
       tc = basepars.lc0_tc; 
       % 2. convolve with the time course
       K_slen = 15;
       n_sti_fr = size(stimpars.x,2);
       tmp = fastconv(stimpars.x(crop_idx_j,:),repmat(tc,K_slen^2,1),K_slen^2,n_sti_fr,basepars.padval);
       c_reg = zeros(15);
       c_reg(6:10,6:10) = 1;
       c_reg = c_reg(:);
       tmp = tmp(c_reg==1,:);
       % 3. zero-mean
       tmp = tmp - repmat(mean(tmp,2),1,n_sti_fr);
       % 4. 
       tmp_p = zeros(size(tmp));
       tmp_n = zeros(size(tmp));
       tmp_p(tmp>0) = tmp(tmp>0);
       tmp_n(tmp<0) = tmp(tmp<0);
       tmp_p = tmp_p.^2;
       tmp_n = tmp_n.^2;
       for k = 1:25
          kx(:,j) = basepars.resc(k)*p(2+k)*tmp_p(k,:)'+kx(:,j);
          kx(:,j) = basepars.resc(k+25)*p(2+25+k)*tmp_n(k,:)'+kx(:,j);
       end
       
    elseif strcmp(basepars.meta_filt_mode,'lc0')
       Kfilt = p(2)*basepars.K;
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       %-- additional non-linear terms
       tc = basepars.lc0_tc;
       % 2. convolve with the time course
       K_slen = 15;
       n_sti_fr = size(stimpars.x,2);
       tmp = fastconv(stimpars.x(crop_idx_j,:),repmat(tc,K_slen^2,1),K_slen^2,n_sti_fr,basepars.padval);
       c_reg = zeros(15);
       c_reg(6:10,6:10) = 1;
       c_reg = c_reg(:);
       tmp = tmp(c_reg==1,:);
       % 3. zero-mean and square
       tmp = tmp - repmat(mean(tmp,2),1,n_sti_fr);
       tmp = tmp.^2;
       for k = 1:25
          kx(:,j) = p(2+k)*tmp(k,:)'+kx(:,j);
       end

    elseif strcmp(basepars.meta_filt_mode,'gbwrk2c')
       Kfilt = p(2)*basepars.K;
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       
       %-- contrast dependencies
       zc = stimpars.dt*std(stimpars.x(crop_idx_j,:))';
       m_zc = mean(zc);
       for k = 1:30
          zc0 = m_zc*ones(1,size(zc,1));
          zc0(k:end) = zc(1:end-k+1);
          kx(:,j) = p(2+k)*zc0'+kx(:,j);
       end

    elseif strcmp(basepars.meta_filt_mode,'gbwrk2mc')
       Kfilt = p(2)*basepars.K;
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       
       %-- mean and contrast dependencies
       zm = stimpars.dt*mean(stimpars.x(crop_idx_j,:))';
       m_zm = mean(zm);
       zc = stimpars.dt*std(stimpars.x(crop_idx_j,:))';
       m_zc = mean(zc);
       for k = 1:30
          zm0 = m_zm*ones(1,size(zm,1));
          zm0(k:end) = zm(1:end-k+1);
          kx(:,j) = p(2+k)*zm0'+kx(:,j);
          
          zc0 = m_zc*ones(1,size(zc,1));
          zc0(k:end) = zc(1:end-k+1);
          kx(:,j) = p(2+30+k)*zc0'+kx(:,j);
       end
       
    elseif strcmp(basepars.meta_filt_mode,'ir6mc')
       for rnk = 1:6
          if rnk == 1
             Kfilt = p(2)*basepars.K{1};
          else
             Kfilt = Kfilt + p(rnk+1)*basepars.K{rnk};
          end
       end
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       
       %-- mean and contrast dependencies
       zm = stimpars.dt*mean(stimpars.x(crop_idx_j,:))';
       m_zm = mean(zm);
       zc = stimpars.dt*std(stimpars.x(crop_idx_j,:))';
       m_zc = mean(zc);
       for k = 1:30
          zm0 = m_zm*ones(1,size(zm,1));
          zm0(k:end) = zm(1:end-k+1);
          kx(:,j) = p(1+6+k)*zm0'+kx(:,j);
          
          zc0 = m_zc*ones(1,size(zc,1));
          zc0(k:end) = zc(1:end-k+1);
          kx(:,j) = p(1+6+30+k)*zc0'+kx(:,j);
       end

    elseif strcmp(basepars.meta_filt_mode,'stc2')
       for rnk = 1:6
          if rnk == 1
             Kfilt = p(2)*basepars.K{1};
          else
             Kfilt = Kfilt + p(rnk+1)*basepars.K{rnk};
          end
       end
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
       
       % --- additinal term
       STC2 = basepars.STC2;
       zstc1 = stimpars.dt*sqrt((stimpars.x(crop_idx_j,:)'*STC2(:,1)).^2)';
       m_zstc1 = mean(zstc1);
       zstc2 = stimpars.dt*sqrt((stimpars.x(crop_idx_j,:)'*STC2(:,2)).^2)';
       m_zstc2 = mean(zstc2);
       
       for k = 1:30
          zstc1t = m_zstc1*ones(1,size(zstc1,2));
          zstc1t(k:end) = zstc1(1:end-k+1);
          kx(:,j) = p(1+6+k)*zstc1t'+kx(:,j);
          
          zstc2t = m_zstc2*ones(1,size(zstc2,2));
          zstc2t(k:end) = zstc2(1:end-k+1);
          kx(:,j) = p(1+6+k)*zstc2t'+kx(:,j);
       end
       
    elseif strcmp(basepars.meta_filt_mode,'rk6mc')
       for rnk = 1:6
          Kfilt = Kfilt + p(rnk+1)*basepars.K{rnk};
       end
       kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
    
    else
       
       switch(basepars.filtermode)
          case {'sta','stark1','gbwrk2'} % edoi
             Kfilt = p(2)*basepars.K;
             kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
             
          case 'indrk2' % edoi
             Kfilt = p(2)*basepars.K{1} + p(3)*basepars.K{2};
             kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,15^2,T,basepars.padval),1); % edoi
             
          case {'nonsep','raw'}
             %keyboard; %--ED
             Kfilt = reshape(p(kidx),basepars.n,basepars.Mk);
             kx(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:),Kfilt,basepars.n,T,basepars.padval),1);
             if (nargout>1)
                1;
                Kfiltsq = reshape(p(ksqidx),basepars.n,basepars.Mk);
                %kxsq(:,j) = stimpars.dt .* sum(fastconv(stimpars.x(crop_idx_j,:).^2,Kfiltsq,basepars.n,T,basepars.padval.^2),1);
                kxsq(:,j) = stimpars.dt .* sum(fastconv(sqrt((stimpars.x(crop_idx_j,:)-0.5).^2),Kfiltsq,basepars.n,T),1);
             end
          case {'sep_raw','rk2'}
             
             [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,basepars);
             
             % Comp1
             kx1 = stimpars.dt .* fastconv(reshape(p(s1_idx),1,basepars.n)*stimpars.x(crop_idx_j,:),reshape(p(t1_idx),1,basepars.Mk),1,T,basepars.padval);
             % Comp 2
             kx2 = stimpars.dt .* fastconv(reshape(p(s2_idx),1,basepars.n)*stimpars.x(crop_idx_j,:),reshape(p(t2_idx),1,basepars.Mk),1,T,basepars.padval);
             kx(:,j) = kx1(:) + kx2(:);


             if (nargout>1)
                [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(ksqidx(1)-1,basepars);
                %kx1 = stimpars.dt .* fastconv(reshape(p(s1_idx),1,basepars.n)*stimpars.x(crop_idx_j,:).^2,reshape(p(t1_idx),1,basepars.Mk),1,T,basepars.padval^2);
                %kx2 = stimpars.dt .* fastconv(reshape(p(s2_idx),1,basepars.n)*stimpars.x(crop_idx_j,:).^2,reshape(p(t2_idx),1,basepars.Mk),1,T,basepars.padval^2);
                kx1 = stimpars.dt .* fastconv(reshape(p(s1_idx),1,basepars.n)*sqrt((stimpars.x(crop_idx_j,:)-0.5).^2),reshape(p(t1_idx),1,basepars.Mk),1,T);
                kx2 = stimpars.dt .* fastconv(reshape(p(s2_idx),1,basepars.n)*sqrt((stimpars.x(crop_idx_j,:)-0.5).^2),reshape(p(t2_idx),1,basepars.Mk),1,T);
                
                
                kxsq(:,j) = kx1(:) + kx2(:);
             end
             
          case 'sep_basis'
             
             
             [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,basepars);
             s1 = basepars.kspace_basis*p(s1_idx);
             t1 = basepars.ktime_basis*p(t1_idx);
             s2 = basepars.kspace_basis*p(s2_idx);
             t2 = basepars.ktime_basis*p(t2_idx);
             
             kx1 = stimpars.dt .* fastconv(s1'*stimpars.x(crop_idx_j,:),t1',1,T,basepars.padval);
             kx2 = stimpars.dt .* fastconv(s2'*stimpars.x(crop_idx_j,:),t2',1,T,basepars.padval);
             kx(:,j) = kx1(:) + kx2(:);
             
             if (nargout>1)
                
                
                [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(ksqidx(1)-1,basepars);
                s1 = basepars.kspace_basis*p(s1_idx);
                t1 = basepars.ktime_basis*p(t1_idx);
                s2 = basepars.kspace_basis*p(s2_idx);
                t2 = basepars.ktime_basis*p(t2_idx);
                
                %kx1 = stimpars.dt .* fastconv(s1'*stimpars.x(crop_idx_j,:).^2,t1',1,T,basepars.padval^2);
                %kx2 = stimpars.dt .* fastconv(s2'*stimpars.x(crop_idx_j,:).^2,t2',1,T,basepars.padval^2);
                kx1 = stimpars.dt .* fastconv(s1'*sqrt((stimpars.x(crop_idx_j,:)-0.5).^2),t1',1,T);
                kx2 = stimpars.dt .* fastconv(s2'*sqrt((stimpars.x(crop_idx_j,:)-0.5).^2),t2',1,T);
                
                kxsq(:,j) = kx1(:) + kx2(:);
                
             end
             
          case 'fixfilt'
             % K = p(get_pars_idx(basepars,j,Neff,'k')).*basepars.K{j};
             % kx(:,j) = stimpars.dt .*  sum(fastconv(stimpars.x(crop_idx_j,:),K,basepars.n,T,basepars.padval),1);
             % modified, edoi, 2011-01-04.
             
             %keyboard;
             %[s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,basepars);
             
             % further modified, 2011-01-05.
             [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,basepars.wnfit.basepars);
             p0K = basepars.wnfit.pstar;
             % Comp1
             %kx1 = stimpars.dt .* fastconv(reshape(p(s1_idx),1,basepars.n)*stimpars.x(crop_idx_j,:),reshape(p(t1_idx),1,basepars.Mk),1,T,basepars.padval);
             kx1 = stimpars.dt .* fastconv(reshape(p0K(s1_idx),1,basepars.n)*stimpars.x(crop_idx_j,:),reshape(p0K(t1_idx),1,basepars.Mk),1,T,basepars.padval);
             % Comp 2
             %kx2 = stimpars.dt .* fastconv(reshape(p(s2_idx),1,basepars.n)*stimpars.x(crop_idx_j,:),reshape(p(t2_idx),1,basepars.Mk),1,T,basepars.padval);
             kx2 = stimpars.dt .* fastconv(reshape(p0K(s2_idx),1,basepars.n)*stimpars.x(crop_idx_j,:),reshape(p0K(t2_idx),1,basepars.Mk),1,T,basepars.padval);
             kx(:,j) = kx1(:) + kx2(:);
             
             if exist('fl_init','var')
                kx = 1*kx;
             else
                kx = p(2)*kx;
             end
             
             if (nargout>1)
                [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(ksqidx(1)-1,basepars);
                %kx1 = stimpars.dt .* fastconv(reshape(p(s1_idx),1,basepars.n)*stimpars.x(crop_idx_j,:).^2,reshape(p(t1_idx),1,basepars.Mk),1,T,basepars.padval^2);
                %kx2 = stimpars.dt .* fastconv(reshape(p(s2_idx),1,basepars.n)*stimpars.x(crop_idx_j,:).^2,reshape(p(t2_idx),1,basepars.Mk),1,T,basepars.padval^2);
                kx1 = stimpars.dt .* fastconv(reshape(p(s1_idx),1,basepars.n)*sqrt((stimpars.x(crop_idx_j,:)-0.5).^2),reshape(p(t1_idx),1,basepars.Mk),1,T);
                kx2 = stimpars.dt .* fastconv(reshape(p(s2_idx),1,basepars.n)*sqrt((stimpars.x(crop_idx_j,:)-0.5).^2),reshape(p(t2_idx),1,basepars.Mk),1,T);
                kxsq(:,j) = kx1(:) + kx2(:);
             end
             
          case 'fixrk2' % edoi, 2011-01-06
             [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,basepars.wnfit.basepars);
             p0K = basepars.wnfit.pstar;
             % Comp1
             kx1 = stimpars.dt .* fastconv(reshape(p0K(s1_idx),1,basepars.n)*stimpars.x(crop_idx_j,:),reshape(p0K(t1_idx),1,basepars.Mk),1,T,basepars.padval);
             % Comp 2
             kx2 = stimpars.dt .* fastconv(reshape(p0K(s2_idx),1,basepars.n)*stimpars.x(crop_idx_j,:),reshape(p0K(t2_idx),1,basepars.Mk),1,T,basepars.padval);
             
             %kx(:,j) = kx1(:) + kx2(:);
             if exist('fl_init','var')
                fprintf('note: filtered stimulus is returned in struct, each filter component separately.\n')
                kx{1} = kx1(:);
                kx{2} = kx2(:);
             else
                kx = p(2)*kx1(:) + p(3)*kx2(:);
             end
             
             if (nargout>1)
                [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(ksqidx(1)-1,basepars);
                %kx1 = stimpars.dt .* fastconv(reshape(p(s1_idx),1,basepars.n)*stimpars.x(crop_idx_j,:).^2,reshape(p(t1_idx),1,basepars.Mk),1,T,basepars.padval^2);
                %kx2 = stimpars.dt .* fastconv(reshape(p(s2_idx),1,basepars.n)*stimpars.x(crop_idx_j,:).^2,reshape(p(t2_idx),1,basepars.Mk),1,T,basepars.padval^2);
                kx1 = stimpars.dt .* fastconv(reshape(p(s1_idx),1,basepars.n)*sqrt((stimpars.x(crop_idx_j,:)-0.5).^2),reshape(p(t1_idx),1,basepars.Mk),1,T);
                kx2 = stimpars.dt .* fastconv(reshape(p(s2_idx),1,basepars.n)*sqrt((stimpars.x(crop_idx_j,:)-0.5).^2),reshape(p(t2_idx),1,basepars.Mk),1,T);
                kxsq(:,j) = kx1(:) + kx2(:);
             end
             
       end
       
    end
end
