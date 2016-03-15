function params = MakeStimKernels(params)
%Gonzalo Mena, 3/2016
Difs=params.patternInfo.Difs;
Diags=params.patternInfo.Diags;
listCurrents=params.patternInfo.listCurrents;
Art=params.patternInfo.Art;
patternNo=params.patternInfo.patternNo;
var0=params.patternInfo.var0;

DifsStim{1}=Difs{1};
DifsStim{2}=zeros(size(Difs{3},1));
options=params.global.options;
breakpoints=findBreakStimElecs(listCurrents);
params.patternInfo.breakpoints=breakpoints;

breakpoints{1}=[0 breakpoints{1}' size(Art,1)];

DiagsStim{1}=Diags{1};
DiagsStim{2}=breakpoints{1}';


for k=1:length(breakpoints{1})-1
DifsStim{2}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1))=Difs{3}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1));
if(DiagsStim{2}(k+1)-DiagsStim{2}(k)>1)
DifsStim{2}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1))=DifsStim{2}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1))/max(max(DifsStim{2}(breakpoints{1}(k)+1:breakpoints{1}(k+1),breakpoints{1}(k)+1:breakpoints{1}(k+1))));
end

end

nVar = [3 length(DiagsStim{2})-1];
nVarCum=cumsum([0 nVar]);
x01=params.arrayInfo.x(1:3);
x02=-4*ones(1,length(DiagsStim{2})-1);
x03=25;
f2 = @(Art,x)logDetKronStimElec(squeeze(Art(:,patternNo,:)),[x log(var0)],DifsStim,[1:3+length(DiagsStim{2})],[1 2],DiagsStim,nVar);
g2 = @(x)f2(Art,x);

xStim = fminunc(g2,[x01 x02 x03],options);


types=[1 2];
for k=1:2
[Ker KerD]=evalKernels(DifsStim{k},DiagsStim{k},xStim(nVarCum(k)+1:nVarCum(k+1)),types(k));
[a b]=eig(Ker);
QStim{k}=a';
QtStim{k}=a;
dLStim{k}=diag(b);
KersStim{k}=Ker;
end

params.patternInfo.KersStim=KersStim;
params.patternInfo.dLStim=dLStim;
params.patternInfo.QtStim=QtStim;
params.patternInfo.QStim=QStim;
params.patternInfo.DifsStim=DifsStim;
params.patternInfo.DiagsStim=DiagsStim;
params.patternInfo.xStim = xStim;