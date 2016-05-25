channel=2;
amplitude=15;
delay=2;

Faza0=[0 0 0 1 2];
Faza1=[0 1 0 1 2];
Faza2=[0 1 1 1 3];
Faza3=[0 1 0 1 1];
Faza4=[0 0 0 1 1];
PhaseData=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
SPD=size(PhaseData);

for p=1:length(PhaseData)
    Phase=PhaseData(p,:)
    [Phase(1:3) de2bi(Phase(4)*amplitude,6,'left-msb') de2bi(Phase(5),4,'left-msb')]
end