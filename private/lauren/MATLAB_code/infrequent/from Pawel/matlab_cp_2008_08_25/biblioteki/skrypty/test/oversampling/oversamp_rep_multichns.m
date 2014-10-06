function y=oversamp_rep_multichns(filenames,spikes,detect_param,filter_param,margins,RC_filter,figures);
size_s=size(filenames);
nrofchns=size_s(2);
spikes_s=size(spikes);
nrofspks=spikes_s(2);

%1 Wykres eenergia bledu vs energia spika
figure(figures(1,1));
hold off;
clf;
hold on;
grid on;
for i=1:nrofchns
    signal=importdata(filenames{i})';
    signal=signal-mean(signal);
    y=blad_vs_energy2(signal,detect_param,filter_param,margins);
    figure(figures(1,1));
    a=plot(y(2,:),y(3,:),'b.');
    kolorow10(i)
    set(a,'Color',kolorow10(i));    
    %hold on;
end
legend(filenames);
xlabel('rms of spike [DAQ units]');
ylabel('rms of interpolation error [DAQ units]');