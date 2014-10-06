% Willems method for calculation of R-D using Shannon Strategy (?)
clc
clear all

Q=3; % total number of items
P=2; % Number of items chosen at a time in J . Que is one item at a time.
% elements of Jcue = 1,2..Q 
% elements of J = 1,2 ... (Q)C(P)
% 
Dict = nchoosek([1:Q]',P);

J_card = size(Dict,1);
J_cue_card = Q;
pjcue=(1/J_cue_card)*ones(J_cue_card,1);
pj_jcue=zeros(J_card,J_cue_card); % J,Jcue
g_true=zeros(J_card,J_cue_card); % input is J, J_cue

p0=0.4;

warning('off','all');

for ji=1:J_card
    for j_cuei=1:J_cue_card
        g_true(ji,j_cuei)=element_in_Dict(ji,j_cuei,Dict);
        pj_jcue(ji,j_cuei) = g_true(ji,j_cuei);
    end
   
end

for j_cuei=1:J_cue_card
 pj_jcue(:,j_cuei) = pj_jcue(:,j_cuei)+0.5; % Arbit!
val = sum(pj_jcue(:,j_cuei))  ;
pj_jcue(:,j_cuei)=pj_jcue(:,j_cuei)/val;
clear val
end

Z_card=2;

pj = pj_jcue*pjcue;

pjcue_j=zeros(J_card,J_card);
for i=1:J_card  % J
    for j=1:J_cue_card% Jcue
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
%%

s_last=0.1;
s_step=0.1;
s_exp=s_last;
resolution=0.01;
s_exp_new=s_exp;
DD_last=0;
iter_ov=0;
while s_exp<2
    s_exp=s_exp_new 
    s_last
    s=-1*(10^s_exp)
   iter_ov=iter_ov+1;
    log_err=[];
    T_card = Z_card^J_cue_card; % Z_card= 3 is number of outputs
    V=log(ones(T_card,J_card)*(1/T_card));  % equivalent to V(t|x)
    
   
    C_sum=zeros(J_card,T_card);
    UB=1;
    LB=-1;
    iter=0;
    while(iter<100000 && UB-LB>0.001) % Some problem with bound convergence right now!
        iter=iter+1;
        C=zeros(J_card,T_card);
        
        % V_Pj_jcue = log(exp(V)*pj_jcue); % typical log-sum-exp problem!
        
        dummy=zeros(T_card,J_cue_card);
        for dumjcue_i=1:J_cue_card
            for dumti=1:T_card
                m=max(V(dumti,:));
                dsumm=0;
                for dumji=1:J_card
                    dsumm = dsumm + exp(V(dumti,dumji)-m)*pj_jcue(dumji,dumjcue_i);
                    
                end
                dummy(dumti,dumjcue_i)=m+log(dsumm);
            end
        end
        if(~exist('V_Pj_jcue'))
        V_Pj_jcue = dummy;
        else
             V_Pj_jcue=  V_Pj_jcue+log(exp(log(2) + dummy -V_Pj_jcue)-1);
        % Dykstra
        end
        
        
        
        for ji=1:J_card
            for ti=1:T_card
                
                summ=0;
                
                
                for j_cuei=1:J_cue_card
                    % if( V(ti,ji)~=0 ||V_Pj_jcue(ti,j_cuei)~=0)
                    summ = summ + pjcue_j(j_cuei,ji)*((V(ti,ji) - V_Pj_jcue(ti,j_cuei)) - s*dist_fcn_2(ji,j_cuei,ti,Z_card,Dict,P,Q)); % doubtful about distortion function!
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
        %V=V_new;
        V=  V +log(exp(log(2) + V_new -V)-1);
        % Dykstra
        exp_V=exp(V);
        C_sum=C_sum+C;
        
        LB = min(C_sum'/iter)*pj; % LB
        % UB
        DD = 0;
        II = 0;
       % Vp_t_jcue = log(exp(V)*pj_jcue); % log sum exp method!! 
      
        Vp_t_jcue  = V_Pj_jcue;
        
        
        
        for ji=1:J_card
            for j_cuei=1:J_cue_card
                for ti=1:T_card
                    DD =DD + pjcue_j(j_cuei,ji)*pj(ji)*exp(V(ti,ji))*dist_fcn_2(ji,j_cuei,ti,Z_card,Dict,P,Q);
                    %wrong II = II+pjcue(j_cuei)*(my_entr(pj_jcue(ji,j_cuei)) + my_entr(V(ti,ji))*pj_jcue(ji,j_cuei));
                    % Use better formulae ?
                    II = II +pjcue(j_cuei)* (my_entr_exp((Vp_t_jcue(ti,j_cuei))) - my_entr_exp((V(ti,ji))))*pj_jcue(ji,j_cuei); % entr of exp .. calculation better ? 
                end
            end
        end
        UB = II - s*DD;
        log_err=[log_err;UB-LB];
        UB-LB;
        DD;
        II;
        exp(V);
        UB
        LB
    end
    DD_log=[DD_log;DD];
    II_log=[II_log;II];
    UB_log=[UB_log;UB];
    LB_log=[LB_log;LB];
    s_log=[s_log;s];
    s_exp_log=[s_exp_log;s_exp];
    p0;
    DD
    DD_last
    iter_ov
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

figure;
plot(DD_log,II_log,'b','Marker','o')
hold on
plot(DD_log,UB_log+s_log.*DD_log,'r','Marker','o')
hold on
plot(DD_log,LB_log+s_log.*DD_log,'g','Marker','o')