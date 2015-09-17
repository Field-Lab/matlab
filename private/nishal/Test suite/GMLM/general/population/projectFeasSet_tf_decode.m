function Z_kp1 = projectFeasSet_tf_decode(ttf,X_kp1,u_k,rho)

T = size(X_kp1,2);
d = size(X_kp1,1);
ttf = [ttf;zeros(T-length(ttf),1)];
fttf = fft(ttf);

figure; idx=0:T-1 ; plot(2*idx/T,abs(fttf));xlabel('freq (*pi) radians');ylabel('amplitude');

% low pass filtering of reconstructed signal
validF = abs(fttf)> max(abs(fttf))*0.05;
freqs = [0:T-1]*2*pi/T;
validFreq = freqs(validF);

E_F = zeros(T,length(validFreq));

for ifreq =1:length(validFreq)
E_F(:,ifreq) = exp(sqrt(-1)*[0:T-1]'*validFreq(ifreq));
E_F(:,ifreq) = E_F(:,ifreq)/sqrt(T);
end
E_F = gpuArray(E_F);
Proj = E_F*E_F';

Z_kp1= 0*X_kp1;
for idim=1:d
x_t = (X_kp1(idim,:) + u_k(idim,:))'; 
Z_kp1(idim,:)= gather(real( E_F*(E_F'*x_t)))';

end

end