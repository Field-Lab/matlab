clear
sigma = load('output/m63_sigma.mat');
sigma = sigma.sigma;
dists = load('output/m63_dists.mat');
dists = dists.dists;
r_mean = load('output/m63_r_mean.mat');
r_mean = r_mean.r_mean;
amp = load('output/m63_amplitude.mat');
amp = amp.amplitude;
offset = load('output/m63_offset.mat');
offset = offset.offset;

%figure(1)
%plot(reshape(dists,1,[]),reshape(sigma,1,[]),'.')

dists(sigma(:)==0) = [];
r_mean(sigma(:)==0) = [];
amp(sigma(:)==0) = [];
offset(sigma(:)==0) = [];
sigma(sigma(:)==0) = [];

dists(abs(sigma(:))>500) = [];
r_mean(abs(sigma(:))>500) = [];
amp(abs(sigma(:))>500) = [];
offset(abs(sigma(:))>500) = [];
sigma(abs(sigma(:))>500) = [];

figure(2)
plot(squeeze(dists),squeeze(sigma),'.')

figure(3)
plot(squeeze(dists),squeeze(amp),'.')

figure(4)
plot(squeeze(dists),squeeze(offset),'.')

figure(5)
plot(squeeze(dists),squeeze(r_mean),'.')