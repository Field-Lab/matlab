function [kx kxsq] = filterstimulus_train2AH(p,Basepars,Stimpars,fl_init)

%%% as of 10-29  atleast rk2 is clean!!! 

%%% Same as filterstimulus_trainAH except got rid of the unnecesarry
%%% Stimpars.dt   term !!! 


%fprintf('\n\n HEY!! filterstimulus_train.m called \n\n')
% FM-DEPENDENT (edoi, 2012-01)



kx = zeros(size(Stimpars.movie_ROI,2),1);
if (nargout>1)
    kxsq = zeros(size(kx));
end
%T= size(Stimpars.movie_ROI,2);
[spacePixels stimFrames] = size(Stimpars.movie_ROI);

    
if (iscell(Basepars.crop_idx))
	crop_idx = Basepars.crop_idx;
else
	crop_idx = Basepars.crop_idx(:,1);
end
    
kidx = Basepars.paramind.L;
if (nargout>1)
    ksqidx = get_pars_idxAH(Basepars,1,1,'ksq');
end
    
    
    
       
switch(Basepars.k_filtermode)
	case {'STA','stark1','gbwrk2'} % edoi
        Kfilt = Basepars.STA;
        kx(:,1) = 1 .* sum(fastconvAH(Stimpars.movie_ROI(crop_idx,:),Kfilt,Basepars.k_spacepixels,stimFrames,Basepars.padval),1); % edoi
	case 'indrk2' % edoi
        Kfilt = p(2)*Basepars.K{1} + p(3)*Basepars.K{2};
        kx(:,1) = 1 .* sum(fastconvAH(Stimpars.movie_ROI(crop_idx,:),Kfilt,Basepars.k_spacepixels,stimFrames,Basepars.padval),1); % edoi
	case {'nonsep','raw'}
    %keyboard; %--ED
        Kfilt = reshape(p(kidx),Basepars.k_spacepixels,Basepars.k_stimframes);
        kx(:,1) = 1 .* sum(fastconvAH(Stimpars.movie_ROI(crop_idx,:),Kfilt,Basepars.k_spacepixels,stimFrames,Basepars.padval),1);
        if (nargout>1)
            1;
            Kfiltsq = reshape(p(ksqidx),Basepars.k_spacepixels,Basepars.k_stimframes);
                %kxsq(:,1) = 1 .* sum(fastconvAH(Stimpars.movie_ROI(crop_idx,:).^2,Kfiltsq,Basepars.k_spacepixels,T,Basepars.padval.^2),1);
            kxsq(:,1) = 1 .* sum(fastconvAH(sqrt((Stimpars.movie_ROI(crop_idx,:)-0.5).^2),Kfiltsq,Basepars.k_spacepixels,stimFrames),1);
        end
    case {'fixedSP'};
        kx = fastconvAH( Stimpars.movie_ROI(crop_idx,:) , (p(kidx)'), 1, stimFrames, Basepars.padval);
    case {'OnOff_hardrect_fixedSP_STA'};
        t1_idx = Basepars.paramind.TIME1 ;
        t2_idx = Basepars.paramind.TIME2 ;       
        kx1 = fastconvAH( Stimpars.movie_ROI(1,:) , (p(t1_idx)'), 1, stimFrames, Basepars.padval);    
        kx2 = fastconvAH( Stimpars.movie_ROI(2,:) , (p(t2_idx)'), 1, stimFrames, Basepars.padval);     
        kx = kx1 + kx2;
     case {'rk2'}            
        s1_idx = Basepars.paramind.SPACE1 ;
        t1_idx = Basepars.paramind.TIME1 ;
        s2_idx = Basepars.paramind.SPACE2 ;
        t2_idx = Basepars.paramind.TIME2 ;
       
        kx1 = fastconvAH( (p(s1_idx)')*Stimpars.movie_ROI(crop_idx,:) , (p(t1_idx)'), 1, stimFrames, Basepars.padval);
        kx2 = fastconvAH( (p(s2_idx)')*Stimpars.movie_ROI(crop_idx,:) , (p(t2_idx)'), 1, stimFrames, Basepars.padval);
        kx  = kx1 + kx2;

        %kx1 = 1 .* fastconvAH(reshape(p(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:),reshape(p(t1_idx),1,Basepars.k_stimframes),1,stimFrames,Basepars.padval);
        %kx2 = 1 .* fastconvAH(reshape(p(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:),reshape(p(t2_idx),1,Basepars.k_stimframes),1,stimFrames,Basepars.padval);
        


        if (nargout>1)
            [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(ksqidx(1)-1,Basepars);
            %kx1 = 1 .* fastconvAH(reshape(p(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:).^2,reshape(p(t1_idx),1,Basepars.k_stimframes),1,T,Basepars.padval^2);
            %kx2 = 1 .* fastconvAH(reshape(p(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:).^2,reshape(p(t2_idx),1,Basepars.k_stimframes),1,T,Basepars.padval^2);
            kx1 = 1 .* fastconvAH(reshape(p(s1_idx),1,Basepars.k_spacepixels)*sqrt((Stimpars.movie_ROI(crop_idx,:)-0.5).^2),reshape(p(t1_idx),1,Basepars.k_stimframes),1,stimFrames);
            kx2 = 1 .* fastconvAH(reshape(p(s2_idx),1,Basepars.k_spacepixels)*sqrt((Stimpars.movie_ROI(crop_idx,:)-0.5).^2),reshape(p(t2_idx),1,Basepars.k_stimframes),1,stimFrames);
            kxsq(:,1) = kx1(:) + kx2(:);
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
             
             kx1 = 1 .* fastconvAH(s1'*Stimpars.movie_ROI(crop_idx,:),t1',1,stimFrames,Basepars.padval);
             kx2 = 1 .* fastconvAH(s2'*Stimpars.movie_ROI(crop_idx,:),t2',1,stimFrames,Basepars.padval);
             kx(:,1) = kx1(:) + kx2(:);
             
             if (nargout>1)
                
                
                [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(ksqidx(1)-1,Basepars);
                s1 = Basepars.kspace_basis*p(s1_idx);
                t1 = Basepars.ktime_basis*p(t1_idx);
                s2 = Basepars.kspace_basis*p(s2_idx);
                t2 = Basepars.ktime_basis*p(t2_idx);
                
                %kx1 = 1 .* fastconvAH(s1'*Stimpars.movie_ROI(crop_idx,:).^2,t1',1,T,Basepars.padval^2);
                %kx2 = 1 .* fastconvAH(s2'*Stimpars.movie_ROI(crop_idx,:).^2,t2',1,T,Basepars.padval^2);
                kx1 = 1 .* fastconvAH(s1'*sqrt((Stimpars.movie_ROI(crop_idx,:)-0.5).^2),t1',1,stimFrames);
                kx2 = 1 .* fastconvAH(s2'*sqrt((Stimpars.movie_ROI(crop_idx,:)-0.5).^2),t2',1,stimFrames);
                
                kxsq(:,1) = kx1(:) + kx2(:);
                
             end
             
          case 'fixfilt'
             % K = p(get_pars_idx(Basepars,1,1,'k')).*Basepars.K{1};
             % kx(:,1) = 1 .*  sum(fastconvAH(Stimpars.movie_ROI(crop_idx,:),K,Basepars.k_spacepixels,T,Basepars.padval),1);
             % modified, edoi, 2011-01-04.
             
             %keyboard;
             %[s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,Basepars);
             
             % further modified, 2011-01-05.
             [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,Basepars.wnfit.Basepars);
             p0K = Basepars.wnfit.pstar;
             % Comp1
             %kx1 = 1 .* fastconvAH(reshape(p(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:),reshape(p(t1_idx),1,Basepars.k_stimframes),1,T,Basepars.padval);
             kx1 = 1 .* fastconvAH(reshape(p0K(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:),reshape(p0K(t1_idx),1,Basepars.k_stimframes),1,stimFrames,Basepars.padval);
             % Comp 2
             %kx2 = 1 .* fastconvAH(reshape(p(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:),reshape(p(t2_idx),1,Basepars.k_stimframes),1,T,Basepars.padval);
             kx2 = 1 .* fastconvAH(reshape(p0K(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:),reshape(p0K(t2_idx),1,Basepars.k_stimframes),1,stimFrames,Basepars.padval);
             kx(:,1) = kx1(:) + kx2(:);
             
             if exist('fl_init','var')
                kx = 1*kx;
             else
                kx = p(2)*kx;
             end
             
             if (nargout>1)
                [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(ksqidx(1)-1,Basepars);
                %kx1 = 1 .* fastconvAH(reshape(p(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:).^2,reshape(p(t1_idx),1,Basepars.k_stimframes),1,T,Basepars.padval^2);
                %kx2 = 1 .* fastconvAH(reshape(p(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:).^2,reshape(p(t2_idx),1,Basepars.k_stimframes),1,T,Basepars.padval^2);
                kx1 = 1 .* fastconvAH(reshape(p(s1_idx),1,Basepars.k_spacepixels)*sqrt((Stimpars.movie_ROI(crop_idx,:)-0.5).^2),reshape(p(t1_idx),1,Basepars.k_stimframes),1,stimFrames);
                kx2 = 1 .* fastconvAH(reshape(p(s2_idx),1,Basepars.k_spacepixels)*sqrt((Stimpars.movie_ROI(crop_idx,:)-0.5).^2),reshape(p(t2_idx),1,Basepars.k_stimframes),1,stimFrames);
                kxsq(:,1) = kx1(:) + kx2(:);
             end
             
          case 'fixrk2' % edoi, 2011-01-06
             [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(kidx(1)-1,Basepars.wnfit.Basepars);
             p0K = Basepars.wnfit.pstar;
             % Comp1
             kx1 = 1 .* fastconvAH(reshape(p0K(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:),reshape(p0K(t1_idx),1,Basepars.k_stimframes),1,stimFrames,Basepars.padval);
             % Comp 2
             kx2 = 1 .* fastconvAH(reshape(p0K(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:),reshape(p0K(t2_idx),1,Basepars.k_stimframes),1,stimFrames,Basepars.padval);
             
             %kx(:,1) = kx1(:) + kx2(:);
             if exist('fl_init','var')
                fprintf('note: filtered stimulus is returned in struct, each filter component separately.\n')
                kx{1} = kx1(:);
                kx{2} = kx2(:);
             else
                kx = p(2)*kx1(:) + p(3)*kx2(:);
             end
             
             if (nargout>1)
                [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(ksqidx(1)-1,Basepars);
                %kx1 = 1 .* fastconvAH(reshape(p(s1_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:).^2,reshape(p(t1_idx),1,Basepars.k_stimframes),1,T,Basepars.padval^2);
                %kx2 = 1 .* fastconvAH(reshape(p(s2_idx),1,Basepars.k_spacepixels)*Stimpars.movie_ROI(crop_idx,:).^2,reshape(p(t2_idx),1,Basepars.k_stimframes),1,T,Basepars.padval^2);
                kx1 = 1 .* fastconvAH(reshape(p(s1_idx),1,Basepars.k_spacepixels)*sqrt((Stimpars.movie_ROI(crop_idx,:)-0.5).^2),reshape(p(t1_idx),1,Basepars.k_stimframes),1,stimFrames);
                kx2 = 1 .* fastconvAH(reshape(p(s2_idx),1,Basepars.k_spacepixels)*sqrt((Stimpars.movie_ROI(crop_idx,:)-0.5).^2),reshape(p(t2_idx),1,Basepars.k_stimframes),1,stimFrames);
                kxsq(:,1) = kx1(:) + kx2(:);
             end
             
       end
       
end

    

%{
%%% DELETED CODE LINES THAT I DON'T USE
%if (strcmp(Basepars.k_filtermode,'fixrk2') && exist('fl_init','var'))
%   kx = cell(2,1);
%else
if (~isfield(Basepars,'k_meta_filt_mode') && strcmp(Basepars.k_filtermode,'stark1'))
   Basepars.k_meta_filt_mode = 'starkn';
end
%}