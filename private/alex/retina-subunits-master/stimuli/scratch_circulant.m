%%
[U S V] = svd((I-P)'*Q'*Q*(I-P));

%%
scatter(fit1.locs_c(:,1),fit1.locs_c(:,2),50,U(:,3),'filled');
set(gca,'YDir','reverse');
set(gca,'CLim',[-0.2 0.2]);
axis off; colormap(cjet);

%%

x1 = randn(32,1);

%%
x1 = mvnrnd(zeros(32,1),M'*M,1)';
x1 = zeros(32,1)/5;
%%
%C = (I-P)'*Q'*Q*(I-P);
%C = (I-P_full)'*Q_full'*Q_full*(I-P_full);
%C = Q_full'*Q_full-M_full'*M_full;
%C = 2.1*Q'*Q - M'*M;
%T = cholcov(C);
x2 = mvnrnd(x1(:,1),C,100);
clf
plot(x2','LineWidth',1);
hold on
plot(x1,'LineWidth',2);
%%
n = 100;
vr_q = mean((repmat(Q*x1,[1,n]) - Q*x2').^2,1);
vr_m = mean((repmat(M*x1,[1,n]) - M*x2').^2,1);
mean(vr_q./vr_m)


%%
for i=1:50
scatter(fit1.locs_c(:,1),fit1.locs_c(:,2),50,x2(i,:),'filled');
set(gca,'YDir','reverse');
set(gca,'CLim',[-2 2]);
axis off; colormap(cjet);
pause(0.1);
end


%%
sz = 32;
[pyr,pind] = buildGpyr(zeros(1,sz));
projection = zeros(size(pyr,1),sz);
for pos=1:sz
  [pyr,pind] = buildGpyr(mkImpulse([1 sz], [1 pos]),'auto','binom5','circular');
  projection(:,pos) = pyr;
end
clf; showIm(projection);

%%
for iband=1:size(pind,1)
    inds = pyrBandIndices(pind,iband);
    bandmat{iband} = projection(inds,:);
end

%% 
Q = bandmat{3};
M = bandmat{4};
P = M'*inv(M*M')*M;
I = eye(32);
for i=1:size(M,1)
    M(i,:) = M(i,:) / norm(M(i,:));
end
for i=1:size(Q,1)
    Q(i,:) = Q(i,:) / norm(Q(i,:));
end
%%
M_block = 8;
sz = 32;
M = [ones(1,M_block) zeros(1,sz-M_block)];
for i=2:sz;
    M(i,:) = circshift(M(i-1,:),[0 1]);
end

Q_block = 4;
Q = [ones(1,Q_block) zeros(1,sz-Q_block)];
for i=2:sz;
    Q(i,:) = circshift(Q(i-1,:),[0 1]);
end

Q = Q(1:Q_block:end,:);
M = M(1:M_block:end,:);

%%
P = M'*inv(M*M')*M;
I = eye(32);


%%
x = linspace(-4,4,32);
%M = normpdf(x,0,1.9);
M = M(1,:)/norm(M(1,:));
for i=2:32;
    M(i,:) = circshift(M(i-1,:),[0 1]);
end

%Q = normpdf(x,0,0.25);
Q = Q(1,:)/norm(Q(1,:));
for i=2:32;
    Q(i,:) = circshift(Q(i-1,:),[0 1]);
end

%M = M(1:8:end,:);
%Q = Q(1:4:end,:);
%P = M'*inv(M*M')*M;
I = eye(32);

imagesc(Q'*Q);

%%
K = [1 zeros(1,31)];
for i=2:32;
    K(i,:) = circshift(K(i-1,:),[0 1]);
end
%%
[U Dq V] = svd(Q'*Q);
[U Dm V] = svd(M'*M);

M = M / sqrt((diag(Dm(1,1))));
Q = Q / sqrt((diag(Dq(1,1))));
%%
C = Q'*Q - M'*M;
[U Dq V] = svd(Q'*Q);
[U Dm V] = svd(M'*M);
[U Dc V] = svd(C);
C = C / (sum(diag(Dc)));
