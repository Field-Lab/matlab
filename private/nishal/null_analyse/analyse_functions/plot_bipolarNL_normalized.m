function plot_bipolarNL_normalized(W,mm,nl_fcn,cols,su_inp_given,su_inp_give)

if(su_inp_given)
    su_inp = su_inp_give;
else
su_inp = W'*mm;
end

range = [-2:0.1:2];
%figure;
for isu=1:size(W,2)
mn = mean(su_inp(isu,:));
sd = sqrt(var(su_inp(isu,:)));
[a,b] = hist(su_inp(isu,:),100); 
bnl = nl_fcn(range*sd+mn);
[hAx,hLine1,hLine2] =plotyy(range,bnl,(b-mn)/sd,a);
hLine1.Color = cols(1);
hLine2.Color = cols(2);
%plot(range,bnl,cols(1));
hold on;

end

end