%Pomiary liniowosci wykonane 8 grudnia 2003. Sygna³ z generatora Agilent 33120A podawany bezpo¶rednio na plytke, gdzie znajduje sie dzielnik napiecia na rezystorach 10kOhm/50Ohm. Wykoanana kalibracja (pomiar z duzym oscyloskopem) wykazala poziom sygnalu 440uV P-P za dzielnikiem przy nastawie sygnalu na generatorze na 1V (amplituda). Pomiary dla sygnalu sinusoidalnego 420 Hz oraz prostokatnego 420Hz i 140Hz. Odczyt sygnalu na wyjsciu za pomoca karty ADC i dalej FFT (program w LabView external_square.vi). Pomiar dla czestosci podstawowej i dwoch harmonicznych (trzykrotnosc i pieciokrotnosc czest. podstawowej). 
%Sygnal byl podawany na wszystkie kanaly rownoczesnie, pomiar na wyjsciu kanalu 15. Plytka numer 1 (z zamontowanymi dodatkowymi potencjometrami).

input=[0.1:0.1:3.5]*440e-6;
input_short=input(1,1:30);

%a) pomiar z sinusem 420 Hz
sin420(1,:)=[22 44 66 88 110 132 154 176 198 220 241 263 284 305 326 347 368 388 409 428 448 467 486 505 523 540 557 574 589 604 619 632 645 657 669]/1000;
sin420=sin420./input*2;

%b) pomiar z prostokatem 420 Hz, wspolczynniki dla 420, 1260, 2100 Hz
prost420(1,:)=[22 45 66 89 111 132 155 176 198 220 241 262 283 304 325 345 364 383 402 420 437 454 469 484 498 513 526 537 547 557]./input_short*2;
prost420(2,:)=[21 41 61 82 102 123 143 163 183 203 223 242 263 282 301 320 338 357 375 392 409 426 442 457 471 487 499 512 523 534]./input_short*2;
prost420(3,:)=[18 35 53 71 88 105 123 141 158 175 192 209 227 244 262 279 295 312 329 345 361 376 391 406 421 437 450 463 475 486]./input_short*2;
prost420=prost420/1000;

%c) pomiar z prostokatem 140Hz, wspolczynniki dla 140, 420, 700Hz
prost140(1,:)=[19 39 58 77 96 115 135 154 172 191 209 227 246 263 280 297 314 329 345 360 374 387 399 411 421 433 443 453 460 468]./input_short*2;
prost140(2,:)=[22 44 66 88 110 132 153 174 195 216 236 255 275 293 311 328 343 358 372 385 396 406 415 423 430 436 442 446 450 454]./input_short*2;
prost140(3,:)=[21 44 66 88 109 130 152 173 194 214 233 252 272 290 307 324 339 354 367 380 391 401 410 418 425 432 438 443 449 454]./input_short*2;
prost140=prost140/1000;

figure(2);
a=plot(input,sin420,input_short,prost420(1,:),input_short,prost420(2,:),input_short,prost420(3,:),input_short,prost140(1,:),input_short,prost140(2,:),input_short,prost140(3,:));
set(a(2),'Marker','d');
set(a(3),'Marker','*');
set(a(4),'Marker','.');

set(a(5),'Marker','d');
set(a(6),'Marker','*');
set(a(7),'Marker','.');
set(a(5),'LineWidth',2.5);
set(a(6),'LineWidth',2.5);
set(a(7),'LineWidth',2.5);

xlabel('input P-P');
ylabel('gain');

legend('420 (sinus)','420 (prost.420)','1260 (prost.420)','2100 (prost. 420)','140 (prost. 140)','420 (prost.140)','700 (prost. 140)');

grid on;
