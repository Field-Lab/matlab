function el_final=NS512_CorrectElectrodesNumbers(electrodes,ArrayClockwiseRotation);

if ArrayClockwiseRotation==270
    el_final=electrodes+128-floor((electrodes+127)/512)*512;
else
    el_final=electrodes;
end
    