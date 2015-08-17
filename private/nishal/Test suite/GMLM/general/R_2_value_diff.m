function [R2_log1,R2_log2] = R_2_value_diff(x1,y1,y2)    

df = abs(y1-y2);
thr = prctile(df,90)
tms = df>thr;

%  
y1 = y1(tms);
y2=y2(tms);
x1=x1(tms);

% R2 value method 2
       n=length(y1);
       r = (n*x1'*y1 - sum(x1)*sum(y1))/(sqrt(n*sum(x1.^2) - sum(x1)^2) * sqrt(n*sum(y1.^2) - sum(y1)^2));
       R2_log1 = r^2;

       n=length(y2);
       r = (n*x1'*y2 - sum(x1)*sum(y2))/(sqrt(n*sum(x1.^2) - sum(x1)^2) * sqrt(n*sum(y2.^2) - sum(y2)^2));
       R2_log2 = r^2;
end