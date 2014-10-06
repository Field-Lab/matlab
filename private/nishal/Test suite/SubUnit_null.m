%% 2 Subunits, 7 cones
nCones=7;

subUnits=[1,1,1,1,0,0,0,0;
          0,0,0,0,1,1,1,1];
      
x=[+1,+1,+1,+1,+1,-5/3,-5/3,-5/3]'/2;

P=perms(x);

subUnitInput = subUnits*(P');

scatter(subUnitInput(1,:)',subUnitInput(2,:)')

%% 4 subunits , 7 cones


nCones=7;
subUnits = [1,1,1,1,0,0,0,0;
            0,0,0,0,1,0,0,0;
            0,0,0,0,0,1,1,0;
            0,0,0,0,0,0,0,1];
      
%x=-[+1,+1,+1,+1,+1,-5/3,-5/3,-5/3]'/2;
x=[+5/3,5/3,5/3,-1,-1,-1,-1,-1]'/2;

P=unique(perms(x),'rows');

subUnitInput = subUnits*(P');
figure;
subplot(3,2,1);
scatter(subUnitInput(1,:)',subUnitInput(2,:)');

subplot(3,2,2);
scatter(subUnitInput(2,:)',subUnitInput(3,:)');


subplot(3,2,3);
scatter(subUnitInput(3,:)',subUnitInput(4,:)');


subplot(3,2,4);
scatter(subUnitInput(1,:)',subUnitInput(3,:)');

subplot(3,2,5);
scatter(subUnitInput(1,:)',subUnitInput(4,:)');

subplot(3,2,6);
scatter(subUnitInput(2,:)',subUnitInput(4,:)');

x=subUnitInput;


%% 
P=double(dec2bin([0:255]))-double('0');
P=(P-0.5)*2;


subUnitInput = subUnits*(P');
figure;
subplot(3,2,1);
scatter(subUnitInput(1,:)',subUnitInput(2,:)');

subplot(3,2,2);
scatter(subUnitInput(2,:)',subUnitInput(3,:)');


subplot(3,2,3);
scatter(subUnitInput(3,:)',subUnitInput(4,:)');


subplot(3,2,4);
scatter(subUnitInput(1,:)',subUnitInput(3,:)');

subplot(3,2,5);
scatter(subUnitInput(1,:)',subUnitInput(4,:)');

subplot(3,2,6);
scatter(subUnitInput(2,:)',subUnitInput(4,:)');
y=subUnitInput;

%%
figure;
subplot(2,1,1);
hist(x(:)',20);
title('Null')

subplot(2,1,2);
hist(y(:)',20);
title('IID white')

%% Cases where atleast one sub-unit OFF
a=(x<0);
b=sum(a,1);
sum(b~=0)/length(b)

a=(y<0);
b=sum(a,1);
sum(b~=0)/length(b)

%% Distribution of number of sub-units ON
figure;
subplot(2,1,1);
a=(x<0);
b=sum(a,1);
[X,N]=hist(b(:));
bar(N,X/sum(X))
ylim([0,1]);
%xlim([0,4])
title('Null')

subplot(2,1,2);
a=(y<0);
b=sum(a,1);
[X,N]=hist(b(:));
bar(N,X/sum(X))
ylim([0,1]);
%xlim([0,4])
title('IID White')

%%

