
proj_log=[];

indx=[59:1:movieLen];
binnedResponses = spkCondColl(condNoMov).spksColl;
icnt=0;
for iTrial=1:nTrials

binnedResponsesTrial=binnedResponses(iTrial,:);

framesValid = indx(binnedResponsesTrial(indx)>0);

for iframe=framesValid
   
    a = movie(:,:,iframe:-1:iframe-Filtlen+1);
    a=a(logical(repmat(Mask,[1,1,Filtlen])));
    xx=a;
%  mask!
proj1 = ((xx'*component1));
proj2 = ((xx'*component2));
icnt=icnt+1;
proj_log(1,icnt) = proj1;
proj_log(2,icnt) = proj2;

end
end


plot(proj_log(1,:),proj_log(2,:),'.','MarkerSize',4.9)
hold on;
plot([0,0],[-2,4],'g');
hold on
plot([-2,4],[0,0],'g');
hold on;
plot(mean(proj_log(1,:)),mean(proj_log(2,:)),'r*');
axis equal