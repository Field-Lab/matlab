% R-D function for different tasks .. just calculate!

% Discrimination task
J_card = 2;
U_card = 3;
pj=(1/J_card)*ones(J_card,1);
pj_jcue=zeros(J_card,J_card); % J,Jcue
g_true=zeros(J_card,J_card); % input is J, J_cue
pjcue = pj; % all equiprobable
p0=0.4;

warning('off','all');

for i=1:J_card
    for j=1:J_card
        g_true(i,j)=sign(i-j);
        pj_jcue(i,j) = (1-p0)*double(i==j) + (p0/(J_card-1))*double(i~=j);
    end
end






D=p0;
%D_log = [0:0.05:p0];
% for D=[0:0.05:p0]
obj_fin_log=[];
isolved=[];
for i=1:3^(U_card*J_card)
    fz = zeros(U_card,J_card); % input is U,J_cue
    
    num=dec2base(i,3,U_card*J_card);
    icount=0;
    for ui = 1:U_card
        for ji=1:J_card
            icount=icount+1;
            fz(ui,ji)=str2double(num(icount)) - 1;
            
        end
    end
    
    
    
    
    % initialize in a constraint feasible way
    
    
    cvx_begin quiet
    variables pu_j(U_card,J_card)
    subject to
    d=0;
    for ji=1:J_card
        for j_cuei=1:J_card
            for ui=1:U_card
                d=d+ pj(j)*pj_jcue(j_cuei,j)*pu_j(ui,ji)*double(g_true(ji,j_cuei)~=fz(ui,j_cuei));
                
            end
        end
    end
    d<=D
    obj = 0;
    pu_j>=0;
    ones(U_card,1)'*pu_j == ones(J_card,1)';
    minimize 0
    cvx_end
    
    if(strcmp(cvx_status,'Solved'))
        pu_jcue = pu_j*pj_jcue; %  initialize!
        
        obj_log=[];
        for iterations_run=1:10
            
            cvx_begin quiet
            variables pu_j(U_card,J_card)
            subject to
            d=0;
            for ji=1:J_card
                for j_cuei=1:J_card
                    for ui=1:U_card
                        d=d+ pj(j)*pj_jcue(j_cuei,j)*pu_j(ui,ji)*double(g_true(ji,j_cuei)~=fz(ui,j_cuei));
                        
                    end
                end
            end
            
            d<=D
            obj = 0;
            
            pu_j>=0;
            ones(U_card,1)'*pu_j == ones(J_card,1)';
            
            
            % Hu_jcue = sum(entr(pu_jcue)*pjcue);
            
            %Hu_j = sum(entr(pu_j)*pj);
            cue_sum=0;
            for j_cuei=1:J_card
                ssum = 0;
                
                for ji=1:J_card
                    ssum_comp=-sum(entr(pu_j(:,ji))) - sum(pu_j(:,ji).*log(pu_jcue(:,j_cuei)));
                    ssum=ssum + pj_jcue(ji,j_cuei)*ssum_comp;
                end
                cue_sum = cue_sum + pjcue(j_cuei)*ssum;
            end
            
            minimize cue_sum%(Hu_jcue - Hu_j)
            cvx_end
            
          
            
            pu_jcue = pu_j*pj_jcue;
            obj_log=[obj_log;cvx_optval];
        end
        plot(obj_log);
        obj_fin_log = [obj_fin_log;obj_log(end)]; 
        isolved = [isolved;i];
          i
          cue_sum
    end

end
cue_sum
%end