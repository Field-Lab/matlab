function [mov_new_ft_re_iftime, mov_new_ft_im_iftime]= proj_null_fourier(sta_mov_ft_iftime,n_cell,mov_ft_re_iftime,mov_ft_im_iftime)

var64=(size(sta_mov_ft_iftime,1)+1)/2;

   A=zeros(2*n_cell,numel(sta_mov_ft_iftime(:,:,1,1))*2);
   s_re=zeros((var64*2-1)*(32*2-1),2*n_cell);
   s_im=zeros((var64*2-1)*(32*2-1),2*n_cell);
    for icell=1:n_cell*2
   
        if(mod(icell,2)==1)
            if(var64==64)
            s_re(:,icell) = reshape(real(sta_mov_ft_iftime(:,:,1,floor((icell-1)/2)+1)),[(64*2-1)*(32*2-1),1]); %s_re{icell}=s_re{icell}(:);
            s_im(:,icell) = reshape(imag(sta_mov_ft_iftime(:,:,1,floor((icell-1)/2)+1)),[(64*2-1)*(32*2-1),1]); %s_im{icell}=s_im{icell}(:);
            end
            
            if(var64==32)
            s_re(:,icell) = reshape(real(sta_mov_ft_iftime(:,:,1,floor((icell-1)/2)+1)),[(32*2-1)*(32*2-1),1]); %s_re{icell}=s_re{icell}(:);
            s_im(:,icell) = reshape(imag(sta_mov_ft_iftime(:,:,1,floor((icell-1)/2)+1)),[(32*2-1)*(32*2-1),1]); %s_im{icell}=s_im{icell}(:);
            end
            A(icell,:) =[s_re(:,icell)',-s_im(:,icell)'];
        else
            A(icell,:)=[s_im(:,icell-1)',s_re(:,icell-1)'];
        end
    end
    m_real = mov_ft_re_iftime(:,:,1);
    m_im = mov_ft_im_iftime(:,:,1);
    m=[m_real(:);m_im(:)];
    
   % null_A = (eye(size(A,2)) - A'*((A*A')\A));
    m=m - A'*((A*A')\(A*m));
    
    mov_new_ft_re_iftime = reshape(m(1:length(m)/2),size(m_real));
    mov_new_ft_im_iftime = reshape(m(length(m)/2+1:end),size(m_im));
   % norm(A*m)
end