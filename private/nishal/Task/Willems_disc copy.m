% Willems method for calculation of R-D using Shannon Strategy (?)
clc
clear all

J_card = 2;
U_card = 3;
pjcue=(1/J_card)*ones(J_card,1);
pj_jcue=zeros(J_card,J_card); % J,Jcue
g_true=zeros(J_card,J_card); % input is J, J_cue

p0=0.2;

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
log_err=[];
T_card = Z_card^J_card; % Z_card= 3 is number of outputs
V=ones(T_card,J_card)*(1/T_card);  % equivalent to V(t|x)
s=-2.5;
C_sum=zeros(J_card,T_card);
UB=1;
LB=-1;
iter=0;
while(UB-LB>0.001)
    iter=iter+1;
    C=zeros(J_card,T_card);
    
    for ji=1:J_card
        for ti=1:T_card
            
            summ=0;
            V_Pj_jcue = V*pj_jcue;
            for j_cuei=1:J_card
               % if( V(ti,ji)~=0 ||V_Pj_jcue(ti,j_cuei)~=0)
                summ = summ + pjcue_j(j_cuei,ji)*(log(V(ti,ji) / (V_Pj_jcue(ti,j_cuei))) - s*dist_fcn(ji,j_cuei,ti,J_card,Z_card)); % doubtful about distortion function! 
               % end
            end
            C(ji,ti) = summ;
        end
    end
    V_Pj_jcue
    C
    
    V_new = 0*V;
    for ji=1:J_card
        V_denj = exp(-C)*V;
        V_denj=V_denj(ji,ji);
        
        for ti=1:T_card
            
        V_new(ti,ji) = V(ti,ji)*exp(-C(ji,ti)) / V_denj;
        end
    end
    V=V_new
    C_sum=C_sum+C;
    
LB = min(C_sum'/iter)*pj; % LB
% UB
DD = 0;
II = 0;
Vp_t_jcue = V*pj_jcue;
for ji=1:J_card
    for j_cuei=1:J_card
        for ti=1:T_card
       DD =DD + pjcue_j(j_cuei,ji)*pj(ji)*V(ti,ji)*dist_fcn(ji,j_cuei,ti,J_card,Z_card); 
       %wrong II = II+pjcue(j_cuei)*(my_entr(pj_jcue(ji,j_cuei)) + my_entr(V(ti,ji))*pj_jcue(ji,j_cuei));
        II = II +pjcue(j_cuei)* (my_entr(Vp_t_jcue(ti,j_cuei)) - my_entr(V(ti,ji)))*pj_jcue(ji,j_cuei);
        end
    end
end
UB = II - s*DD;
log_err=[log_err;UB-LB];
UB-LB
DD
II
V
UB
LB
end

p0
DD     
II
