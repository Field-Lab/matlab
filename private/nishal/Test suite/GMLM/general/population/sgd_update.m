
function z = sgd_update(validSpectrum,rho,lambda,init_pt)

z=init_pt;
for iter=1:10000
[~,~,z] = power_spectrum_match(z,validSpectrum,rho,lambda,init_pt,0.1/iter);
end

end