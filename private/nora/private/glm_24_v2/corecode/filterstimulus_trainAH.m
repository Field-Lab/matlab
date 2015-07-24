function [kx kxsq] = filterstimulus_trainAH(p,Basepars,Stimpars,fl_init)

%fprintf('\n\n HEY!! filterstimulus_train.m called \n\n')
% FM-DEPENDENT (edoi, 2012-01)

Neff = 1;
if (strcmp(Basepars.k_filtermode,'fixrk2') && exist('fl_init','var'))
   kx = cell(2,1);
else
   kx = zeros(size(Stimpars.movie_ROI,2),Neff);
end
if (nargout>1)
    kxsq = zeros(size(kx));
end

T = size(Stimpars.movie_ROI,2);

if (~isfield(Basepars,'k_meta_filt_mode') && strcmp(Basepars.k_filtermode,'stark1'))
   Basepars.k_meta_filt_mode = 'starkn';
end


for j=1:1
    
    if (iscell(Basepars.crop_idx))
        crop_idx_j = Basepars.crop_idx{j};
    else
        crop_idx_j = Basepars.crop_idx(:,j);
    end
    
    kidx = get_pars_idxAH(Basepars,j,Neff,'k');
    if (nargout>1)
        ksqidx = get_pars_idxAH(Basepars,j,Neff,'ksq');
    end
    
    
    
       
    switch(Basepars.k_filtermode)
          case {'sta','stark1','gbwrk2'} % edoi
             Kfilt = p(2)*Basepars.K;
             kx(:,j) = Stimpars.dt .* sum(fastconv(Stimpars.movie_ROI(crop_idx_j,:),Kfilt,Basepars.k_spacepixels,T,Basepars.padval),1); % edoi
          case 'indrk2' % edoi
             Kfilt = p(2)*Basepars.K{1} + p(3)*Basepars.K{2};
             kx(:,j) = Stimpars.dt .* sum(fastconv(Stimpars.movie_ROI(crop_idx_j,:),Kfilt,Basepars.k_spacepixels,T,Basepars.padval),1); % edoi
          case {'nonsep','raw'}
             %keyboard; %--ED
             Kfilt = reshape(p(kidx),Basepars.k_spacepixels,Basepars.k_stimframes);
             kx(:,j) = Stimpars.dt .* sum(fastconv(Stimpars.movie_ROI(crop_idx_j,:),Kfilt,Basepars.k_spacepixels,T,Basepars.padval),1);
             if (nargout>1)
                1;
                Kfiltsq = reshape(p(ksqidx),Basepars.k_spacepixels,Basepars.k_stimframes);
                %kxsq(:,j) = Stimpars.dt .* sum(fastconv(Stimpars.movie_ROI(crop_idx_j,:).^2,Kfiltsq,Basepars.k_spacepixels,T,Basepars.padval.^2),1);
                kxsq(:,j) = Stimpars.dt .* sum(fastconv(sqrt((Stimpars.movie_ROI(crop_idx_j,:)-0.5).^2),Kfiltsq,Basepars.k_spacepixels,T),1);
             end
          case {'sep_raw','rk2'}
              
              spacepixels = Basepars.k_spacepixels;
              stimframes = Basepars.k_stimframes;
              
              offset = kidx(1)-1;
              s1_idx = offset+1:offset+spacepixels;
              offset = offset + spacepixels;
              t1_idx = offset+1:offset+stimframes;
              offset = offset + stimframes;
              s2_idx = offset+1:offset+spacepixels;
              offset = offset + spacepixels;
              t2_idx = offset+1:offset+stimframes;
            %%% FOLLOWING IS REPLACED WITH THE ACTUAL MANUAL CODE WRITTEN
            %%% ABOVE 
            %% [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,Basepars);
             
            
             % Comp1
             kx1 = Stimpars.dt .* fastconv(reshape(p(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:),reshape(p(t1_idx),1,Basepars.k_stimframes),1,T,Basepars.padval);
             % Comp 2
             kx2 = Stimpars.dt .* fastconv(reshape(p(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:),reshape(p(t2_idx),1,Basepars.k_stimframes),1,T,Basepars.padval);
             kx(:,j) = kx1(:) + kx2(:);


             if (nargout>1)
                [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(ksqidx(1)-1,Basepars);
                %kx1 = Stimpars.dt .* fastconv(reshape(p(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:).^2,reshape(p(t1_idx),1,Basepars.k_stimframes),1,T,Basepars.padval^2);
                %kx2 = Stimpars.dt .* fastconv(reshape(p(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:).^2,reshape(p(t2_idx),1,Basepars.k_stimframes),1,T,Basepars.padval^2);
                kx1 = Stimpars.dt .* fastconv(reshape(p(s1_idx),1,Basepars.k_spacepixels)*sqrt((Stimpars.movie_ROI(crop_idx_j,:)-0.5).^2),reshape(p(t1_idx),1,Basepars.k_stimframes),1,T);
                kx2 = Stimpars.dt .* fastconv(reshape(p(s2_idx),1,Basepars.k_spacepixels)*sqrt((Stimpars.movie_ROI(crop_idx_j,:)-0.5).^2),reshape(p(t2_idx),1,Basepars.k_stimframes),1,T);
                
                
                kxsq(:,j) = kx1(:) + kx2(:);
             end
             
          case 'sep_basis'
             
             spacepixels = Basepars.nofilters_kspace;
             stimframes  = Basepars.nofilters_ktime;
         %    [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,Basepars);
             
             offset = kidx(1)-1;
             s1_idx = offset+1:offset+spacepixels;
             offset = offset + spacepixels;
             t1_idx = offset+1:offset+stimframes;
             offset = offset + stimframes;
             s2_idx = offset+1:offset+spacepixels;
             offset = offset + spacepixels;
             t2_idx = offset+1:offset+stimframes;
             
             
             s1 = Basepars.kspace_basis*p(s1_idx);
             t1 = Basepars.ktime_basis*p(t1_idx);
             s2 = Basepars.kspace_basis*p(s2_idx);
             t2 = Basepars.ktime_basis*p(t2_idx);
             
             kx1 = Stimpars.dt .* fastconv(s1'*Stimpars.movie_ROI(crop_idx_j,:),t1',1,T,Basepars.padval);
             kx2 = Stimpars.dt .* fastconv(s2'*Stimpars.movie_ROI(crop_idx_j,:),t2',1,T,Basepars.padval);
             kx(:,j) = kx1(:) + kx2(:);
             
             if (nargout>1)
                
                
                [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(ksqidx(1)-1,Basepars);
                s1 = Basepars.kspace_basis*p(s1_idx);
                t1 = Basepars.ktime_basis*p(t1_idx);
                s2 = Basepars.kspace_basis*p(s2_idx);
                t2 = Basepars.ktime_basis*p(t2_idx);
                
                %kx1 = Stimpars.dt .* fastconv(s1'*Stimpars.movie_ROI(crop_idx_j,:).^2,t1',1,T,Basepars.padval^2);
                %kx2 = Stimpars.dt .* fastconv(s2'*Stimpars.movie_ROI(crop_idx_j,:).^2,t2',1,T,Basepars.padval^2);
                kx1 = Stimpars.dt .* fastconv(s1'*sqrt((Stimpars.movie_ROI(crop_idx_j,:)-0.5).^2),t1',1,T);
                kx2 = Stimpars.dt .* fastconv(s2'*sqrt((Stimpars.movie_ROI(crop_idx_j,:)-0.5).^2),t2',1,T);
                
                kxsq(:,j) = kx1(:) + kx2(:);
                
             end
             
          case 'fixfilt'
             % K = p(get_pars_idx(Basepars,j,Neff,'k')).*Basepars.K{j};
             % kx(:,j) = Stimpars.dt .*  sum(fastconv(Stimpars.movie_ROI(crop_idx_j,:),K,Basepars.k_spacepixels,T,Basepars.padval),1);
             % modified, edoi, 2011-01-04.
             
             %keyboard;
             %[s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,Basepars);
             
             % further modified, 2011-01-05.
             [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,Basepars.wnfit.Basepars);
             p0K = Basepars.wnfit.pstar;
             % Comp1
             %kx1 = Stimpars.dt .* fastconv(reshape(p(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:),reshape(p(t1_idx),1,Basepars.k_stimframes),1,T,Basepars.padval);
             kx1 = Stimpars.dt .* fastconv(reshape(p0K(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:),reshape(p0K(t1_idx),1,Basepars.k_stimframes),1,T,Basepars.padval);
             % Comp 2
             %kx2 = Stimpars.dt .* fastconv(reshape(p(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:),reshape(p(t2_idx),1,Basepars.k_stimframes),1,T,Basepars.padval);
             kx2 = Stimpars.dt .* fastconv(reshape(p0K(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:),reshape(p0K(t2_idx),1,Basepars.k_stimframes),1,T,Basepars.padval);
             kx(:,j) = kx1(:) + kx2(:);
             
             if exist('fl_init','var')
                kx = 1*kx;
             else
                kx = p(2)*kx;
             end
             
             if (nargout>1)
                [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(ksqidx(1)-1,Basepars);
                %kx1 = Stimpars.dt .* fastconv(reshape(p(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:).^2,reshape(p(t1_idx),1,Basepars.k_stimframes),1,T,Basepars.padval^2);
                %kx2 = Stimpars.dt .* fastconv(reshape(p(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:).^2,reshape(p(t2_idx),1,Basepars.k_stimframes),1,T,Basepars.padval^2);
                kx1 = Stimpars.dt .* fastconv(reshape(p(s1_idx),1,Basepars.k_spacepixels)*sqrt((Stimpars.movie_ROI(crop_idx_j,:)-0.5).^2),reshape(p(t1_idx),1,Basepars.k_stimframes),1,T);
                kx2 = Stimpars.dt .* fastconv(reshape(p(s2_idx),1,Basepars.k_spacepixels)*sqrt((Stimpars.movie_ROI(crop_idx_j,:)-0.5).^2),reshape(p(t2_idx),1,Basepars.k_stimframes),1,T);
                kxsq(:,j) = kx1(:) + kx2(:);
             end
             
          case 'fixrk2' % edoi, 2011-01-06
             [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,Basepars.wnfit.Basepars);
             p0K = Basepars.wnfit.pstar;
             % Comp1
             kx1 = Stimpars.dt .* fastconv(reshape(p0K(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:),reshape(p0K(t1_idx),1,Basepars.k_stimframes),1,T,Basepars.padval);
             % Comp 2
             kx2 = Stimpars.dt .* fastconv(reshape(p0K(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:),reshape(p0K(t2_idx),1,Basepars.k_stimframes),1,T,Basepars.padval);
             
             %kx(:,j) = kx1(:) + kx2(:);
             if exist('fl_init','var')
                fprintf('note: filtered stimulus is returned in struct, each filter component separately.\n')
                kx{1} = kx1(:);
                kx{2} = kx2(:);
             else
                kx = p(2)*kx1(:) + p(3)*kx2(:);
             end
             
             if (nargout>1)
                [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(ksqidx(1)-1,Basepars);
                %kx1 = Stimpars.dt .* fastconv(reshape(p(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:).^2,reshape(p(t1_idx),1,Basepars.k_stimframes),1,T,Basepars.padval^2);
                %kx2 = Stimpars.dt .* fastconv(reshape(p(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx_j,:).^2,reshape(p(t2_idx),1,Basepars.k_stimframes),1,T,Basepars.padval^2);
                kx1 = Stimpars.dt .* fastconv(reshape(p(s1_idx),1,Basepars.k_spacepixels)*sqrt((Stimpars.movie_ROI(crop_idx_j,:)-0.5).^2),reshape(p(t1_idx),1,Basepars.k_stimframes),1,T);
                kx2 = Stimpars.dt .* fastconv(reshape(p(s2_idx),1,Basepars.k_spacepixels)*sqrt((Stimpars.movie_ROI(crop_idx_j,:)-0.5).^2),reshape(p(t2_idx),1,Basepars.k_stimframes),1,T);
                kxsq(:,j) = kx1(:) + kx2(:);
             end
             
       end
       
    end
end
