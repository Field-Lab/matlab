% Two sub-units
% Assume [1,0] and [0,1] are two sub-units. 

% stim_grid = -2:2
dS=0.01;
stims = -10:dS:10;

sigma=2;
meanStim=0;
probStims =dS*(1/(sqrt(2*pi)*sigma))*(exp(-((stims-meanStim).^2)/(2*sigma^2))); % Gaussian

%probStims = double(stims>0.5 & stims<1.5) + double(stims>-1.5 & stims<-0.5);
probStims=probStims/sum(probStims); % As range truncated bor us! 

Threshold = sum(probStims)
probMap = zeros(length(stims));
respMap = zeros(length(stims));
%hinge = @(x) x.*double(x>0);
%respFunc = @(x,y) exp(hinge(x).^2)+exp(hinge(y).^2);
 respFunc= @(x,y) exp(x) + exp(y);

for i=1:size(probMap,1)
    for j=1:size(probMap,2)
        probMap(i,j)=probStims(j)*probStims(i);
        respMap (i,j) = respFunc(stims(j),stims(i)); %exp(stims(i)) + exp(stims(j));
    end
end

probMap = probMap/sum(probMap(:));


 [X,Y]=meshgrid(stims,stims);
figure('Color','w');
subplot(2,2,1);
contourf(X,Y,probMap/max(probMap(:)),20);
hold on
contour(X,Y,respMap/max(respMap(:)),100);
axis image
hold on
plot([1,0,1/sqrt(2)],[0,1,(1/sqrt(2))],'*')
xlabel('Sub-unit 1');
ylabel('Sub-unit 2');
title('Input stimulus and Iso-response curves');

subplot(2,2,2);
SpikeTriggerdEnsemble = probMap.*respMap;
contourf(X,Y,SpikeTriggerdEnsemble,20);
axis image
STA=zeros(2,1);
STA(1) = sum(stims.*mean(SpikeTriggerdEnsemble,1)/sum(mean(SpikeTriggerdEnsemble,1)));
STA(2) = sum(stims'.*mean(SpikeTriggerdEnsemble,2)/sum(mean(SpikeTriggerdEnsemble,2)));
hold on
plot([STA(1),0,1],[STA(2),1,0],'*')
hold on
plot(0*stims,stims);
plot(stims,0*stims);

xlabel('Sub-unit 1');
ylabel('Sub-unit 2');
title('Spike Triggered Ensemble');

% Response on null line
subplot(2,2,3);
S1 = stims;
S2 = -STA(1)*S1/STA(2);
respNull=respFunc(S1,S2);
plotyy(S1,probStims,S1,respNull);
xlabel('Input to sub-unit 2 on Null axis');
ylabel('Response to Input');
hold on
spikeTriggeredNullResp= probStims.*respNull;
spikeTriggeredNullResp=spikeTriggeredNullResp/sum(spikeTriggeredNullResp);
plot(S1,spikeTriggeredNullResp,'r')
legend('Probability Of Input','Response','Spike Triggered Ensemble')

figure;
[X,Y]=meshgrid(stims,stims);
actualPartition = double(1*X+0*Y>=0*X+1*Y);
plot(stims,stims,'r');
hold on
fittedPartition = double(0.6*X+0.4*Y>=0.4*X+0.6*Y);
%contour(X,Y,fittedPartition);

title('Actual and Fitted Partition');
hold on
contour(X,Y,probMap/max(probMap(:)),20);

Sigma3Boundary = sqrt(X.^2+Y.^2)<3*sigma;
hold on;
contour(X,Y,Sigma3Boundary);
hold on;
plot([0,1],[1,0],'b*');
hold on;
plot([0.6,0.4],[0.4,0.6],'r*');

%% 3D case

dS=0.01;
stims = -1:dS:1;

x1=stims;
x2=stims; 

%respF3D =@(x,y,z)(exp(x)+exp(y)+exp(z));
respF3D =@(x,y,z)((x).*(x>0)+(y).*(y>0)+(z).*(z>0));

probMap = zeros(length(stims));
respMap = zeros(length(stims));
%hinge = @(x) x.*double(x>0);
%respFunc = @(x,y) exp(hinge(x).^2)+exp(hinge(y).^2);
 respFunc= @(x,y) exp(x) + exp(1*y);
 
[X,Y]=meshgrid(stims,stims);

for i=1:size(probMap,1)
    for j=1:size(probMap,2)
        sumP=0;
        for delta = stims
            Pi = probMap(i);
            Pj = probMap(j);
        sumP=sumP+(dS*(1/(sqrt(2*pi)*sigma))*(exp(-((Pi+delta-meanStim).^2)/(2*sigma^2)))) * (dS*(1/(sqrt(2*pi)*sigma))*(exp(-((Pj+delta-meanStim).^2)/(2*sigma^2)))) * (dS*(1/(sqrt(2*pi)*sigma))*(exp(-((-Pi-Pj+delta-meanStim).^2)/(2*sigma^2))));
        end
        
        probMap(i,j)=sumP;
        respMap (i,j) = respF3D(stims(j),stims(i),-(stims(i)+stims(j))); %exp(stims(i)) + exp(stims(j));
        if ((X(i,j)^2+Y(i,j)^2 + (X(i,j)+Y(i,j))^2)>1.1)
          % respMap(i,j)=0;
           probMap(i,j)=0;
        end
        
            if ((X(i,j)^2+Y(i,j)^2 + (X(i,j)+Y(i,j))^2)<0.8)
            probMap(i,j)=0;
            end
    end
end

probMap = probMap / sum(probMap(probMap>0));


figure('Color','w');
subplot(1,3,1);
contourf(X,Y,probMap/max(probMap(:)),20);
hold on
%contour(X,Y,respMap/max(respMap(:)),100);
axis image
hold on
plot([1,-1,1,-1,0,0]/sqrt(2),[-1,1,0,0,1,-1]/sqrt(2),'r*')
hold on;
plot([1,-2,1]/sqrt(6),[1,1,-2]/sqrt(6),'b*');
hold on;
plot([2,-1,-1]/sqrt(6),[-1,2,-1]/sqrt(6),'g*');

xlabel('Sub-unit 1');
ylabel('Sub-unit 2');
title('Input stimulus and Iso-response curves');

subplot(1,3,2);
SpikeTriggerdEnsemble = probMap.*respMap;
contourf(X,Y,SpikeTriggerdEnsemble,20);
axis image
hold on
plot([1,-1,1,-1,0,0]/sqrt(2),[-1,1,0,0,1,-1]/sqrt(2),'r*')
hold on;
plot([1,-2,1]/sqrt(6),[1,1,-2]/sqrt(6),'b*');
hold on;
plot([2,-1,-1]/sqrt(6),[-1,2,-1]/sqrt(6),'g*');
