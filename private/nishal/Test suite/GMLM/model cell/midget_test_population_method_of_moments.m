% Test population 
% 2 midget cells share a sub-unit..

nCones = 14;
nSU = 8;
nCells = 2;

SU_cone_mat = zeros(nSU,nCones);
SU_cone_mat(1,[1,2])=1;
SU_cone_mat(2,[3,4])=1;
SU_cone_mat(3,[5,6])=1;
SU_cone_mat(4,[7])=1;
SU_cone_mat(5,[8,9])=1;
SU_cone_mat(6,[10,11])=1;
SU_cone_mat(7,12)=1;
SU_cone_mat(8,[13,14])=1;
SU_cone_wt = SU_cone_mat;

cell_su_mat = zeros(nCells,nSU);
cell_su_mat(1,[1,2,3,4,5])=1;
cell_su_mat(2,[5,6,7,8])=1;
cell_su_wt = 2*cell_su_mat;

% f =@(x)
 f=@(x)exp(x);
%f = @(x) max(x/2,0);
g=@(x)x;

figure;
subplot(1,2,1);
imagesc(imresize(SU_cone_mat,10,'nearest'));
axis image
set(gca,'xTick',[]);
set(gca,'yTick',[]);
xlabel('Cones');
ylabel('SU');

subplot(1,2,2);
imagesc(imresize(cell_su_wt,10,'nearest'));
axis image
set(gca,'xTick',[]);
set(gca,'yTick',[]);
xlabel('SU');
ylabel('Cells');
%% 
T = 120*120*60;
sig=2/3; %2/5
movie_c = sig*randn(nCones,T);
su_inp = SU_cone_wt*movie_c;
cell_op = g(cell_su_wt*f(su_inp));
dt=1/120;

figure;
[x,n] = hist(su_inp(1,:),30);
plotyy(n,x,n,f(n));


resp = poissrnd(cell_op*dt);
firingRate = sum(resp,2)/(120*60);

%% empirical and true estimates 

icell=1;
moment_emp = zeros(nCones,nCones);
moment = moment_emp;
gamma=2;
for ipix = 1:nCones
    
    for jpix = 1:nCones
    
        if(ipix==jpix)
       
            u=zeros(nCones,1);
            u(1)=gamma; u(2)=-gamma;
            moment_emp(ipix,jpix) =  sum(resp(icell,:)'/T ) * sum(exp(u'*movie_c))/T; 
            moment(ipix,jpix) = sum(2*exp((sum(SU_cone_wt(logical(cell_su_mat(icell,:)'),:).^2,2)*sig^2 )/2),1)*dt * exp((sum(u.^2)*sig^2) / 2); 
      
        else
            u=zeros(nCones,1);
            u(ipix)=gamma; u(jpix)=-gamma;
            moment_emp(ipix,jpix) = exp(u'*movie_c)*resp(icell,:)'/T; 
            
            
            xx=SU_cone_wt(logical(cell_su_mat(icell,:)'),:);
            xx = xx + repmat(u',[size(xx,1),1]);
            moment(ipix,jpix) = sum(2*exp((sum(xx.^2,2)*sig^2 )/2),1) *dt; 
        end 
    end
end

moment_diff = moment - repmat(diag(moment),[1,nCones]); 
moment_emp_diff = moment_emp - repmat(diag(moment_emp),[1,nCones]);

figure;

subplot(1,2,1);
imagesc(imresize(moment_diff,10,'nearest'));
colorbar
axis image
title('True statistic difference');
set(gca,'xTick',[]);
set(gca,'yTick',[]);

subplot(1,2,2);
imagesc(imresize(moment_emp_diff,10,'nearest'));
title('Estimated statistic difference');
colorbar
axis image

set(gca,'xTick',[]);
set(gca,'yTick',[]);



%% Different weights
SU_cone_wt(1,1) = 0.7;
SU_cone_wt(1,2) = 1;
icell=1;

 gamma_expt = cell(3,1); gamma_emp_expt=cell(3,1);
gamma_range = 0.5:0.5:5;
 
 for jjpix=[2,3]
     
 gamma1_eq_log = [];
 gamma2_eq_log = [];
 
 gamma1_emp_eq_log = [];
 gamma2_emp_eq_log = [];
for gamma_norm = gamma_range;
    gamma_norm
    
gamma1_log=[];moment_emp_11_log=[];moment_emp_12_log=[];moment_11_log=[];moment_12_log=[];

for e=0:0.01:gamma_norm^2
    gamma1 = sqrt(e);
   gamma2 = -sqrt(gamma_norm^2-gamma1^2); 

   ipix=1;jpix=1;
  u=zeros(nCones,1);
  u(1)=gamma1; u(2)=gamma2;
  moment_emp_11 =  sum(resp(icell,:)'/T ) * sum(exp(u'*movie_c))/T; 
  moment_11 = sum(2*exp((sum(SU_cone_wt(logical(cell_su_mat(icell,:)'),:).^2,2)*sig^2 )/2),1)*dt * exp((sum(u.^2)*sig^2) / 2); 
      
  ipix=1;jpix=jjpix;
   u=zeros(nCones,1);
   u(ipix)=gamma1; u(jpix)=gamma2;
   moment_emp_12 = exp(u'*movie_c)*resp(icell,:)'/T; 
   
   xx=SU_cone_wt(logical(cell_su_mat(icell,:)'),:);
   xx = xx + repmat(u',[size(xx,1),1]);
   moment_12 = sum(2*exp((sum(xx.^2,2)*sig^2 )/2),1) *dt; 

   gamma1_log = [gamma1_log;gamma1];
   moment_emp_11_log = [moment_emp_11_log;moment_emp_11];
   moment_11_log = [moment_11_log;moment_11];
   moment_emp_12_log = [moment_emp_12_log;moment_emp_12];
   moment_12_log = [moment_12_log;moment_12];
end

 figure;
 plot((gamma1_log.^2)/(gamma_norm^2),moment_12_log);
 hold on;
 plot((gamma1_log.^2)/(gamma_norm^2),moment_11_log)

 title(sprintf('Pixels 1-%d : gamma norm : %f',jpix,gamma_norm));
% hold on;plot((gamma1_log.^2)/(gamma_norm^2),moment_emp_12_log);hold on;plot((gamma1_log.^2)/(gamma_norm^2),moment_emp_11_log)

 
 
 % intersection -moment
 [v,I] = min((moment_11_log-moment_12_log).^2);
 
 gamma1_eq = gamma1_log(I);
 gamma2_eq = -sqrt(gamma_norm^2-gamma1_eq^2); 
 
 gamma1_eq_log = [gamma1_eq_log;gamma1_eq];
 gamma2_eq_log = [gamma2_eq_log;gamma2_eq];
 
 
 
 % intersection -moment empirical
 [v,I] = min((moment_emp_11_log-moment_emp_12_log).^2);
 
 gamma1_eq = gamma1_log(I);
 gamma2_eq = -sqrt(gamma_norm^2-gamma1_eq^2); 
 
 gamma1_emp_eq_log = [gamma1_emp_eq_log;gamma1_eq];
 gamma2_emp_eq_log = [gamma2_emp_eq_log;gamma2_eq];
 
 
end

figure; 
plot(gamma2_eq_log,gamma1_eq_log); hold on;
plot([0,0],[-max(abs(gamma1_eq_log)),max(abs(gamma1_eq_log))],'g');
hold on;
plot([-max(abs(gamma2_eq_log)),max(abs(gamma2_eq_log))],[0,0],'g');
hold on;
plot(gamma2_emp_eq_log,gamma1_emp_eq_log,'--');



gamma_expt{jjpix} = [gamma1_eq_log,gamma2_eq_log];
gamma_emp_expt{jjpix} = [gamma1_emp_eq_log,gamma2_emp_eq_log];
 end
 
 
figure; 
plot(gamma_expt{2}(:,2),gamma_expt{2}(:,1),'r');
hold on;
plot(gamma_expt{3}(:,2),gamma_expt{3}(:,1)),'b';
hold on;
plot(gamma_emp_expt{2}(:,2),gamma_emp_expt{2}(:,1),'--r');
hold on;
plot(gamma_emp_expt{3}(:,2),gamma_emp_expt{3}(:,1),'--b');
hold on;
plot([0,0],[-max(abs(gamma1_eq_log)),max(abs(gamma1_eq_log))],'g');
hold on;
plot([-max(abs(gamma2_eq_log)),max(abs(gamma2_eq_log))],[0,0],'g');


figure; 
plot(gamma_expt{2}(:,2)./gamma_range',gamma_expt{2}(:,1)./gamma_range','--r*');
hold on;
plot(gamma_expt{3}(:,2)./gamma_range',gamma_expt{3}(:,1)./gamma_range','--b.');
hold on;
plot(gamma_emp_expt{2}(:,2)./gamma_range',gamma_emp_expt{2}(:,1)./gamma_range','--m*');
hold on;
plot(gamma_emp_expt{3}(:,2)./gamma_range',gamma_emp_expt{3}(:,1)./gamma_range','--c.');
hold on;
plot([0,0],[-max(abs(gamma1_eq_log)),max(abs(gamma1_eq_log))],'g');
hold on;
plot([-max(abs(gamma2_eq_log)),max(abs(gamma2_eq_log))],[0,0],'g');

%% Different weights -imaginary probe
SU_cone_wt(1,1) = 0.7;
SU_cone_wt(1,2) = 1;
icell=1;

 gamma_expt = cell(3,1); gamma_emp_expt=cell(3,1);
gamma_range = 0.5:0.5:5;
 
 for jjpix=[2,3]
     
 gamma1_eq_log = [];
 gamma2_eq_log = [];
 
 gamma1_emp_eq_log = [];
 gamma2_emp_eq_log = [];
for gamma_norm = gamma_range;
    gamma_norm
    
gamma1_log=[];moment_emp_11_log=[];moment_emp_12_log=[];moment_11_log=[];moment_12_log=[];

for e=0:0.01:gamma_norm^2
    gamma1 = sqrt(e);
   gamma2 = -sqrt(gamma_norm^2-gamma1^2);
   
   
  u=zeros(nCones,1);
  u(1)=gamma1; u(2)=gamma2; u = sqrt(-1)*u;
  moment_emp_11 =  sum(resp(icell,:)'/T ) * sum(exp(u'*movie_c))/T; 
  moment_11 = sum(2*exp((sum(SU_cone_wt(logical(cell_su_mat(icell,:)'),:).^2,2)*sig^2 )/2),1)*dt * exp((sum(u.^2)*sig^2) / 2); 
      
  ipix=1;jpix=jjpix;
   u=zeros(nCones,1);
   u(ipix)=gamma1; u(jpix)=gamma2;u = sqrt(-1)*u;
   moment_emp_12 = exp(u'*movie_c)*resp(icell,:)'/T; 
   
   xx=SU_cone_wt(logical(cell_su_mat(icell,:)'),:);
   xx = xx + repmat(u',[size(xx,1),1]);
   moment_12 = sum(2*exp((sum(xx.^2,2)*sig^2 )/2),1) *dt; 

   gamma1_log = [gamma1_log;gamma1];
   moment_emp_11_log = [moment_emp_11_log;moment_emp_11];
   moment_11_log = [moment_11_log;moment_11];
   moment_emp_12_log = [moment_emp_12_log;moment_emp_12];
   moment_12_log = [moment_12_log;moment_12];
end

 figure;
 plot((gamma1_log.^2)/(gamma_norm^2),real(moment_12_log));
 hold on;
 plot((gamma1_log.^2)/(gamma_norm^2),real(moment_11_log));
 hold on;
 plot((gamma1_log.^2)/(gamma_norm^2),imag(moment_12_log));
 hold on;
 plot((gamma1_log.^2)/(gamma_norm^2),imag(moment_11_log));

 title(sprintf('Pixels 1-%d : gamma norm : %f',jpix,gamma_norm));
% hold on;plot((gamma1_log.^2)/(gamma_norm^2),moment_emp_12_log);hold on;plot((gamma1_log.^2)/(gamma_norm^2),moment_emp_11_log)

 
 
 % intersection -moment
 [v,I] = min(abs(moment_11_log-moment_12_log).^2);
 
 gamma1_eq = gamma1_log(I);
 gamma2_eq = -sqrt(gamma_norm^2-gamma1_eq^2); 
 
 gamma1_eq_log = [gamma1_eq_log;gamma1_eq];
 gamma2_eq_log = [gamma2_eq_log;gamma2_eq];
 
 
 
 % intersection -moment empirical
 [v,I] = min((moment_emp_11_log-moment_emp_12_log).^2);
 
 gamma1_eq = gamma1_log(I);
 gamma2_eq = -sqrt(gamma_norm^2-gamma1_eq^2); 
 
 gamma1_emp_eq_log = [gamma1_emp_eq_log;gamma1_eq];
 gamma2_emp_eq_log = [gamma2_emp_eq_log;gamma2_eq];
 
 
end

figure; 
plot(gamma2_eq_log,gamma1_eq_log); hold on;
plot([0,0],[-max(abs(gamma1_eq_log)),max(abs(gamma1_eq_log))],'g');
hold on;
plot([-max(abs(gamma2_eq_log)),max(abs(gamma2_eq_log))],[0,0],'g');
hold on;
plot(gamma2_emp_eq_log,gamma1_emp_eq_log,'--');



gamma_expt{jjpix} = [gamma1_eq_log,gamma2_eq_log];
gamma_emp_expt{jjpix} = [gamma1_emp_eq_log,gamma2_emp_eq_log];
 end
 
 
figure; 
plot(gamma_expt{2}(:,2),gamma_expt{2}(:,1),'r');
hold on;
plot(gamma_expt{3}(:,2),gamma_expt{3}(:,1)),'b';
hold on;
plot(gamma_emp_expt{2}(:,2),gamma_emp_expt{2}(:,1),'--r');
hold on;
plot(gamma_emp_expt{3}(:,2),gamma_emp_expt{3}(:,1),'--b');
hold on;
plot([0,0],[-max(abs(gamma1_eq_log)),max(abs(gamma1_eq_log))],'g');
hold on;
plot([-max(abs(gamma2_eq_log)),max(abs(gamma2_eq_log))],[0,0],'g');


figure; 
plot(gamma_expt{2}(:,2)./gamma_range',gamma_expt{2}(:,1)./gamma_range','--r*');
hold on;
plot(gamma_expt{3}(:,2)./gamma_range',gamma_expt{3}(:,1)./gamma_range','--b.');
hold on;
plot(gamma_emp_expt{2}(:,2)./gamma_range',gamma_emp_expt{2}(:,1)./gamma_range','--m*');
hold on;
plot(gamma_emp_expt{3}(:,2)./gamma_range',gamma_emp_expt{3}(:,1)./gamma_range','--c.');
hold on;
plot([0,0],[-max(abs(gamma1_eq_log)),max(abs(gamma1_eq_log))],'g');
hold on;
plot([-max(abs(gamma2_eq_log)),max(abs(gamma2_eq_log))],[0,0],'g');