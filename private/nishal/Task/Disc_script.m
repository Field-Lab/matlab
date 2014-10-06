%%
matlabpool
%%
Z_card=3;
J_card_log=cell(4*9,1);
p0_log=cell(4*9,1);
icnt=0;
for J_card=[2,3,4,5]
    for p0=[0.1:0.1:0.9]
        icnt=icnt+1;
        J_card_log{icnt}=J_card;
        p0_log{icnt}=p0;
    %J_card_log=[J_card_log;J_card];
    %p0_log=[p0_log;p0];
    end
end
%%
%matlabpool
parfor i=1:36
    i
Willems_disc_smart_fcn(J_card_log{i},p0_log{i},Z_card)
end

%% 
%Willems_disc_smart_fcn(2,0.2,3)