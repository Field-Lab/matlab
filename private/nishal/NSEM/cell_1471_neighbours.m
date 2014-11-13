
c1=cell(118/2,1);
icnt=1;
for itrial=1:1:118/2
c1{itrial}=CBP_spk1{icnt};
icnt=icnt+2;
end


c2=cell(118/2,1);
icnt=1;
for itrial=1:1:118/2
c2{itrial}=CBP_spk2{icnt};
icnt=icnt+2;
end

c3=cell(118/2,1);
icnt=1;
for itrial=1:1:118/2
c3{itrial}=CBP_spk3{icnt}';
icnt=icnt+2;
end

c4=cell(118/2,1);
icnt=1;
for itrial=1:1:118/2
c4{itrial}=CBP_spk4{icnt}';
icnt=icnt+2;
end

figure;
[x1,y1]=plotSpikeRaster(c1,'PlotType','vertline')
figure;
[x2,y2]=plotSpikeRaster(c2,'PlotType','vertline')
figure;
[x3,y3]=plotSpikeRaster(c3,'PlotType','vertline')
figure;
[x4,y4]=plotSpikeRaster(c4,'PlotType','vertline')

figure;
plot(x1, y1, 'k');
hold on
plot(x2, y2+28, 'r');
hold on
plot(x3, y3+28*2, 'b');
hold on
plot(x4, y4+28*3, 'm');


%%
spkCondColl(4).a=[];
spkCondColl(1).spksColl=CBP_spk1;
spkCondColl(2).spksColl=CBP_spk2;
spkCondColl(3).spksColl=CBP_spk3;
spkCondColl(4).spksColl=CBP_spk4;

% RE FORMAT spkCondColl
nTrials=27;
spkCondCollformat(4).spksColl=[];
for icond=1:4
    spksColl=zeros(7200,nTrials);
    
    for itrial=1:nTrials
    for ispk=1:length(spkCondColl(icond).spksColl{itrial})
       frameNo= floor(double(spkCondColl(icond).spksColl{itrial}(ispk))*120)+1;
       spksColl(frameNo,itrial)=1;
    end
    end
    spkCondCollformat(icond).spksColl=logical(spksColl)';
end


addpath(genpath('../null_analyse/'));

psthBinSize=10;
psthSmoothen=5;
for icond=1:4
[timeLogData{icond},psthData{icond}]=  psth_calc(( spkCondCollformat(icond).spksColl),psthBinSize,'nonoverlap');
psthData{icond}=conv(psthData{icond},(1/psthSmoothen)*ones(psthSmoothen,1),'same');
end


figure;
plot(timeLogData{1}/120,psthData{1},'k');
hold on
plot(timeLogData{2}/120,psthData{2},'r');
hold on
plot(timeLogData{3}/120,psthData{3},'b');
hold on
plot(timeLogData{4}/120,psthData{4},'m');
xlabel('Time (s)')


%% 

% Refractory violations ..

ref_viol=zeros(4,4);

for icell1=1:4
for icell2=1:4
num_viol=0;
tot_spk=0;

for itrial=1:30
for ilen=1:length(spkCondColl(icell1).spksColl{itrial})
    
if(icell1==icell2)
    num_viol=num_viol+double(sum(double(abs(spkCondColl(icell1).spksColl{itrial}(ilen)-spkCondColl(icell2).spksColl{itrial})<=0.002))-1>0);
end
if(icell1~=icell2)
    num_viol=num_viol+double(sum(double(abs(spkCondColl(icell1).spksColl{itrial}(ilen)-spkCondColl(icell2).spksColl{itrial})<=0.002))>0);
end
end

tot_spk=tot_spk+length(spkCondColl(icell1).spksColl{itrial});

end
ref_viol(icell1,icell2)=100*num_viol/tot_spk;
end
end
