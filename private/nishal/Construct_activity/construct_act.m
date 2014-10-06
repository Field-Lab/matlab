

k = ggsim.k; % input stimulus
k=k(end:-1:1);
k_len=length(k);
h=ggsim.ih; % Post spike filter
h_len=length(h);
dt=ggsim.dt;

time_len=20; %s
time_sample_len=time_len/dt;


%% 
%lam_des = (ones(time_sample_len,1).*(3+30*double(mod([1:time_sample_len]',100)>50)));
lam_des=4*ones(time_sample_len,1);
plot(lam_des);
log_lam_des=log(lam_des);

K=zeros(time_sample_len,time_sample_len);
H=zeros(time_sample_len,time_sample_len);
for itime=1:time_sample_len
    ik_len=min(k_len,itime);%min(itime+k_len,time_sample_len)-itime;
K(itime,itime-ik_len+1:itime) = k(1:ik_len)';

    ih_len=min(itime+h_len,time_sample_len)-itime;
H(itime:itime+ih_len-1,itime) = h(1:ih_len);
end

%% Use only 1 time poi 
cvx_begin 
variables x(time_sample_len,1) y(time_sample_len,1)

subject to
norm(y,1)+norm(1-y,1)<=time_sample_len
%x<=2
%x>=-2
minimize (norm(log_lam_des-K*x-H*y-log(4),2) - 0.01* y'*log_lam_des + 0.001*norm(x))
cvx_end

%% Independent time model!

%Replace y by lam_des .. for norm calculation, its wrong some variance, etc
%.. 
cvx_begin 
variables x(time_sample_len,1) y(time_sample_len,1)

subject to
%norm(y,1)+norm(1-y,1)<=time_sample_len
%x<=2
%x>=-2
minimize (norm(log_lam_des-K*x-H*lam_des*dt-log(4),2) - 0.1* dt*lam_des'*log_lam_des )
cvx_end
%%

x = (K'*K)\K'*(log_lam_des-H*lam_des*dt-log(4));

%%

%Replace y by lam_des .. for norm calculation, its wrong some variance, etc
%.. 
Ht = H(:,1:end-1); % Last is column of zeros!
Htinv = sparse((Ht'*Ht)\Ht');

cvx_begin 
variables x(time_sample_len,1) y(time_sample_len,1)

subject to
%norm(y,1)+norm(1-y,1)<=time_sample_len
%x<=2
%x>=-2

q= Htinv*(log_lam_des-K*x-log(4));
minimize ((q'*q -2*q'*lam_des(1:end-1)*dt + sum(lam_des(1:end-1)*dt)) - 0.1* dt*lam_des'*log_lam_des )
cvx_end


%%
y_exp_log=[];
for itrial=1:50

      y_exp=[];
      lin_inp=K*x;
    for itime=1:time_sample_len
      
    log_lam_exp = lin_inp(itime) +H(itime,1:itime-1)*y_exp + log(4);
    if(rand(1)<exp(log_lam_exp)*dt)
    y_exp=[y_exp;1];
    else
    y_exp=[y_exp;0];
    end
    

    end
y_exp_log=[y_exp_log;y_exp'];    
end

mean_spks=mean(y_exp_log,1);
% PSTH ? 
time_bin_len=10;
PSTh=zeros(time_sample_len,1);
for itime=1:time_sample_len-time_bin_len
   PSTH(itime)=sum(mean_spks(itime:itime+time_bin_len))/(dt*time_bin_len); 
end


tms=[];
for itrial=1:50
    for itime=1:size(y_exp_log,2)
        if(y_exp_log(itrial,itime)==1)
        tms=[tms,(itrial-1)*size(y_exp_log,2)+itime]; 
        end
    end
end

rasterplot(tms',50,size(y_exp_log,2));


figure;
plot(PSTH);hold on
plot(lam_des,'r')
legend('PSTH','desired firing rate pattern')

figure;
plot(x(1:1900))
title('Input stimulus')
