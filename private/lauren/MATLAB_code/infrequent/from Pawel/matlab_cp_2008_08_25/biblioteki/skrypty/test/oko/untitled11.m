for i=1:16
subplot(4,4,i);
semilogy(a(1:8192,i));
end

figure(3);
semilogy(a(1:8192,20));

cd /home/pawel/pliki/nauka/neuro64/lipiec2002;