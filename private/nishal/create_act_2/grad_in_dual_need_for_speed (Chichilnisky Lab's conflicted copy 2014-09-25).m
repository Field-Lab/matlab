% Gradient ascent!!


%KKt=filt_mat*filt_mat';
b=-filt_mat*mov_reshape;
alpha=0.1;
beta=0.5;
momentum=0.7;
%%
% Gradient descent is oscillatory .. Do momoentum version??

tic;
nu=zeros(size(filt_mat,1),1);
f_nu=fcn_eval(nu,filt_mat,b);
dx = -(0.5*(filt_mat*(filt_mat'*nu)) + b);


res_log=[];
prev_step=0*nu;
for iter=1:1000
    iter
    
    % Now do line search
    togo=1;
    t=1;
    
while togo==1
    f_ev=fcn_eval(nu+t*dx,filt_mat,b);
    f_comp=f_nu+alpha*t*(-dx'*dx);
    if(f_ev<f_comp)
        togo=0;
        nu_new=nu+t*dx+momentum*prev_step; % I doubt if momentum term need to be added to line search part as well!
        prev_step=nu_new-nu;
        norm(nu_new-nu)
        nu=nu_new;
        f_nu=fcn_eval(nu,filt_mat,b)
        dx = -(0.5*(filt_mat*(filt_mat'*nu)) + b);
    else
        t=beta*t;
    end
end
 


mov_modify=(-0.5)*(filt_mat'*nu)+mov_reshape;
res=norm(filt_mat*mov_modify)



res_log=[res_log;res];
end
toc;


mov_modify=(-0.5)*(filt_mat'*nu)+mov_reshape;
res=norm(filt_mat*mov_modify)