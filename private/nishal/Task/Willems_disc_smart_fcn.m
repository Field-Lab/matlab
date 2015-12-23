% Willems method for calculation of R-D using Shannon Strategy (?)
function []=Willems_disc_smart_fcn(J_card,p0,Z_card)
%J_card = 2;

pjcue=(1/J_card)*ones(J_card,1);
pj_jcue=zeros(J_card,J_card); % J,Jcue
g_true=zeros(J_card,J_card); % input is J, J_cue

%p0=0.2;

warning('off','all');

for i=1:J_card
    for j=1:J_card
        g_true(i,j)=sign(i-j);
        pj_jcue(i,j) = (1-p0)*double(i==j) + (p0/(J_card-1))*double(i~=j);
    end
end
Z_card=3;

pj = pj_jcue*pjcue;

pjcue_j=zeros(J_card,J_card);
for i=1:J_card  % J
    for j=1:J_card % Jcue
        pjcue_j(j,i) = pj_jcue(i,j)*pjcue(j)/pj(i);
    end
end
%%

%s=-3;
DD_log=[];
II_log=[];
UB_log=[];
LB_log=[];
s_log=[];
s_exp_log=[];

s_last=0.0;
s_step=0.1;
s_exp=s_last;
resolution=0.01;
s_exp_new=s_exp;
DD_last=0;
iter_ov=0;
while s_exp<2
    s_exp=s_exp_new; 
    s_last;
    s=-1*(10^s_exp);
   iter_ov=iter_ov+1;
    log_err=[];
    T_card = Z_card^J_card; % Z_card= 3 is number of outputs
    V=log(ones(T_card,J_card)*(1/T_card));  % equivalent to V(t|x)
    
    C_sum=zeros(J_card,T_card);
    UB=1;
    LB=-1;
    iter=0;
    while(iter<1000)%(UB-LB>0.001) % Some problem with bound convergence right now!
        iter=iter+1;
        C=zeros(J_card,T_card);
        
        % V_Pj_jcue = log(exp(V)*pj_jcue); % typical log-sum-exp problem!
        
        dummy=zeros(T_card,J_card);
        for dumjcue_i=1:J_card
            for dumti=1:T_card
                m=max(V(dumti,:));
                dsumm=0;
                for dumji=1:J_card
                    dsumm = dsumm + exp(V(dumti,dumji)-m)*pj_jcue(dumji,dumjcue_i);
                    
                end
                dummy(dumti,dumjcue_i)=m+log(dsumm);
            end
        end
        
        V_Pj_jcue = dummy;
        
        
        
        
        for ji=1:J_card
            for ti=1:T_card
                
                summ=0;
                
                
                for j_cuei=1:J_card
                    % if( V(ti,ji)~=0 ||V_Pj_jcue(ti,j_cuei)~=0)
                    summ = summ + pjcue_j(j_cuei,ji)*((V(ti,ji) - V_Pj_jcue(ti,j_cuei)) - s*dist_fcn(ji,j_cuei,ti,J_card,Z_card)); % doubtful about distortion function!
                    % end
                end
                C(ji,ti) = summ;
            end
        end
        V_Pj_jcue;
        exp_V_Pj_jcue = exp(V_Pj_jcue);
        C;
        
        V_new = ones(size(V));
        for ji=1:J_card
            V_denj = log(exp(-C)*exp(V));
            V_denj=V_denj(ji,ji);
            
            for ti=1:T_card
                
                V_new(ti,ji) = V(ti,ji)-C(ji,ti)-V_denj;
            end
        end
        V=V_new;
        exp_V=exp(V);
        C_sum=C_sum+C;
        
        LB = min(C_sum'/iter)*pj; % LB
        % UB
        DD = 0;
        II = 0;
        Vp_t_jcue = log(exp(V)*pj_jcue);
        for ji=1:J_card
            for j_cuei=1:J_card
                for ti=1:T_card
                    DD =DD + pjcue_j(j_cuei,ji)*pj(ji)*exp(V(ti,ji))*dist_fcn(ji,j_cuei,ti,J_card,Z_card);
                    %wrong II = II+pjcue(j_cuei)*(my_entr(pj_jcue(ji,j_cuei)) + my_entr(V(ti,ji))*pj_jcue(ji,j_cuei));
                    % Use better formulae ?
                    II = II +pjcue(j_cuei)* (my_entr(exp(Vp_t_jcue(ti,j_cuei))) - my_entr(exp(V(ti,ji))))*pj_jcue(ji,j_cuei);
                end
            end
        end
        UB = II - s*DD;
        log_err=[log_err;UB-LB];
        UB-LB;
        DD;
        II;
        exp(V);
        UB;
        LB;
    end
   DD_log=[DD_log;DD];
    II_log=[II_log;II];
    UB_log=[UB_log;UB];
    LB_log=[LB_log;LB];
    s_log=[s_log;s];
    s_exp_log=[s_exp_log;s_exp];
    p0
    DD
    DD_last
    iter_ov
    J_card
    II;
    UB;
    LB;
    if(iter_ov>1)
    if(abs(DD-DD_last)>resolution)
    s_exp_new = (s_exp+s_last)/2;
    else
        while(DD~=min(DD_log)&&~isempty(DD_log((DD-DD_log)<resolution & DD_log<DD)))
        D_stuff= (DD_log.*double((DD-DD_log)<resolution).*double(DD_log<DD));
        D_stuff(D_stuff==0)=Inf;    
        [DD,indxx] = min(D_stuff);
        s_exp_new = s_exp_log(indxx);
        end
            s_last=s_exp_new;
        s_exp_new = s_exp_new+s_step;
        DD_last=DD;
    end
    else 
        s_last=s_exp_new;
        s_exp_new = s_exp_new+s_step;
        DD_last=DD;
    end
    
    
    
end


[s_log,indx]=sort(s_log);
DD_log = DD_log(indx);
UB_log = UB_log(indx);
LB_log = LB_log(indx);
II_log = II_log(indx);

f=figure;
plot(DD_log,II_log,'b','Marker','o')
hold on
plot(DD_log,UB_log+s_log.*DD_log,'r','Marker','o')
hold on
plot(DD_log,LB_log+s_log.*DD_log,'g','Marker','o')
title(sprintf('J %d, p0 %0.02f , Z_card %d',J_card,p0,Z_card));

 print(f,sprintf('/Volumes/Analysis/nishal/Task/J_%d_p0_%0.02f_Z_%d',J_card,p0,Z_card),'-dpng');
save(sprintf('/Volumes/Analysis/nishal/Task/J_%d_p0_%0.02f_Z_%d.mat',J_card,p0,Z_card));
end
