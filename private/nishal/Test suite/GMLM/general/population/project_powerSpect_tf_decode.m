function Z_kp1 = project_powerSpect_tf_decode(ttf,X_kp1,u_k,rho,lambda,scale_spec)

T = size(X_kp1,2);
d = size(X_kp1,1);
ttf = [ttf;zeros(T-length(ttf),1)];
fttf = fft(ttf);

figure; idx=0:T-1 ; plot(2*idx/T,abs(fttf));xlabel('freq (*pi) radians');ylabel('amplitude');

% low pass filtering of reconstructed signal
validF = abs(fttf)> max(abs(fttf))*0;
freqs = [0:T-1]*2*pi/T;
validFreq = freqs(validF);
validSpectrum = (abs(fttf(validF)).^2)*scale_spec; % *scale_spec

optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter-detailed',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','off');

Z_kp1= 0*X_kp1;

tic;
parfor idim=1:d
%  x_init = (X_kp1(idim,:) + u_k(idim,:))'; 
% Z_kp1(idim,:)= gather(real( E_F*(E_F'*x_t)))';
% Z_kp1(idim,:) = fminunc(@(z)power_spectrum_match(z,validSpectrum,rho,lambda,(X_kp1(idim,:) + u_k(idim,:))'),randn(T,1),optim_struct)';

% SGD
Z_kp1(idim,:) = sgd_update(validSpectrum,rho,lambda,(X_kp1(idim,:) + u_k(idim,:))');

end
toc;

end