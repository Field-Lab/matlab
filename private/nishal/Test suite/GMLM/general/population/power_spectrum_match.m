function [fval,grad,update] = power_spectrum_match(z,validSpectrum,rho,lam,x_init,etak)

diff = abs(fft(z)).^2 - validSpectrum;
fval=0;
%fval = (lam/2)*norm(diff)^2 + (rho/2)*norm(z-x_init)^2;

grad = lam*real(2*ifft((diff.*fft(z)))) - rho*(x_init-z); 

update = z - etak*grad;
end 