function w=rgb(x);
s=2^24;

r=floor((s-x)/65536);
g=floor((s-x-r*65536)/256);
b=s-x-r*65536-g*256;

w=[r g b];