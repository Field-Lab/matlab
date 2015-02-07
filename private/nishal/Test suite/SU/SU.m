% Two sub-units
% Assume [1,0] and [0,1] are two sub-units. 

% stim_grid = -2:2
dS=0.01;
stims = -5:dS:5;

sigma=2;
meanStim=0;
probStims =dS*(1/(sqrt(2*pi)*sigma))*(exp(-((stims-meanStim).^2)/(2*sigma^2))); % Gaussian

%probStims = double(stims>0.5 & stims<1.5) + double(stims>-1.5 & stims<-0.5);
probStims=probStims/sum(probStims); % As range truncated bor us! 


probMap = zeros(length(stims));
respMap = zeros(length(stims));
%hinge = @(x) x.*double(x>0);
%respFunc = @(x,y) exp(hinge(x).^2)+exp(hinge(y).^2);
 respFunc= @(x,y) exp(x) + exp(1.1*y);

for i=1:size(probMap,1)
    for j=1:size(probMap,2)
        probMap(i,j)=probStims(i)*probStims(j);
        respMap (i,j) = respFunc(stims(i),stims(j)); %exp(stims(i)) + exp(stims(j));
    end
end
 [X,Y]=meshgrid(stims,stims);
figure;
subplot(1,3,1);
contourf(X,Y,probMap/max(probMap(:)),20);
hold on
contour(X,Y,respMap/max(respMap(:)),20);
axis image
hold on
plot([1,0,1/sqrt(2)],[0,1,(1/sqrt(2))],'*')

subplot(1,3,2);
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

% Response on null line
subplot(1,3,3);
S1 = stims;
S2 = -STA(1)*S1/STA(2);
respNull=respFunc(S1,S2);
plotyy(S1,probStims,S1,respNull);
xlabel('Input to sub-unit 1 on Null axis');
ylabel('Response to Input');
hold on
spikeTriggeredNullResp= probStims.*respNull;
spikeTriggeredNullResp=spikeTriggeredNullResp/sum(spikeTriggeredNullResp);
plot(S1,spikeTriggeredNullResp,'r')
legend('Probability Of Input','Response','Spike Triggered Ensemble')

