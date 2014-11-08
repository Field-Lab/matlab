tic;
A=zeros(2*n_cell,numel(sta_mov_ft{1}(:,:,iftime))*2);
for iftime=movie_len:-1:1
   %iftime
    
   
   s_re=cell(2*n_cell,1);
   s_im=cell(2*n_cell,1);
    for icell=1:n_cell*2
   
        if(mod(icell,2)==1)
            s_re{icell} = real(sta_mov_ft{floor((icell-1)/2)+1}(:,:,iftime)); s_re{icell}=s_re{icell}(:);
            s_im{icell} = imag(sta_mov_ft{floor((icell-1)/2)+1}(:,:,iftime)); s_im{icell}=s_im{icell}(:);
            A(icell,:) =[s_re{icell}',-s_im{icell}'];
        else
            A(icell,:)=[s_im{icell-1}',s_re{icell-1}'];
        end
    end
    m_real = mov_ft_re(:,:,iftime);
    m_im = mov_ft_im(:,:,iftime);
    m=[m_real(:);m_im(:)];
    
   % null_A = (eye(size(A,2)) - A'*((A*A')\A));
    m=m - A'*((A*A')\(A*m));
    
    mov_new_ft_re(:,:,iftime) = reshape(m(1:length(m)/2),size(m_real));
    mov_new_ft_im(:,:,iftime) = reshape(m(length(m)/2+1:end),size(m_im));
   % norm(A*m)
    
end
toc;
