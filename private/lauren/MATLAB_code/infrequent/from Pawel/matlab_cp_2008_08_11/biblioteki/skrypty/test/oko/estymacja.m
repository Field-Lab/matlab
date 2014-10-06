n=4096;
M=64;
N=M*n;

a=rand(1,N)-0.5;

%fN=fft(a);
[fN,wN,covN]=gwsz_f(a,20000,1);
subplot(4,1,1);
plot(fN,wN);

fn=zeros(1,n);

for i=1:M
   s=a(1,(1+(i-1)*n):i*n);
   %fi=fft(s);
   [fi,wi,convi]=gwsz_f(s,20000,1);
   fn=fn+wi;
end
fn=fn/M;

subplot(4,1,2);
plot(wi);
subplot(4,1,3);
plot(fn);

mean(mean(wN))
mean(mean(wi))
mean(mean(fn))

[f,w]=spdf(a,n,n,20000,1);
subplot(4,1,4);
plot(f,w);
mean(mean(fn))
mean(mean(w))
std(fn)
std(w)