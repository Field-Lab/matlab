% Gradient ascent!!


%KKt=filt_mat*filt_mat';

%b=-filt_mat*mov_reshape;
function []=grad_in_dual_need_for_speed(stas,mov,movie_time,n_cell,filt_dim1,filt_dim2,filt_len,mov_number)

b= -Ax(stas,mov,movie_time,n_cell);
alpha=0.1;
beta=0.5;
momentum=0.7;
%%
% Gradient descent is oscillatory .. Do momoentum version??

tic;
%nu=zeros(size(filt_mat,1),1);
nu=zeros(movie_time,n_cell);
f_nu=fcn_eval_new(nu,stas,b)
dx = -(0.5*(Ax(stas,Atx(stas,nu,filt_dim1,filt_dim2,filt_len,movie_time,n_cell),movie_time,n_cell)) + b);


res_log=[];
prev_step=0*nu;
for iter=1:1000
    iter
    
    % Now do line search
    togo=1;
    t=1;
    
while togo==1
    f_ev=fcn_eval_new(nu+t*dx,stas,b);%fcn_eval(nu+t*dx,filt_mat,b);
    f_comp=f_nu+alpha*t*(-sum(dx(:).^2));
    if(f_ev<f_comp)
        togo=0;
        nu_new=nu+t*dx+momentum*prev_step; % I doubt if momentum term need to be added to line search part as well!
        prev_step=nu_new-nu;
        norm(nu_new(:)-nu(:));
        nu=nu_new;
       f_nu=fcn_eval_new(nu,stas,b);
dx = -(0.5*(Ax(stas,Atx(stas,nu,filt_dim1,filt_dim2,filt_len,movie_time,n_cell),movie_time,n_cell)) + b);

    else
        t=beta*t;
    end
end
 


%mov_modify_new=(-0.5)*(Atx(stas,nu,filt_dim1,filt_dim2,filt_len,movie_time,n_cell))+mov;
%mat=Ax(stas,mov_modify_new,movie_time,n_cell);
%res=norm(mat(:))


%res_log=[res_log;res];
end
toc;


mov_modify_new=(-0.5)*(Atx(stas,nu,filt_dim1,filt_dim2,filt_len,movie_time,n_cell))+mov;
%mat=Ax(stas,mov_modify_new,movie_time,n_cell);
save(sprintf('/Volumes/Analysis/nishal/movies/mov_%d',mov_number),'mov_modify_new');
end
