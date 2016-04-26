function response = sig_NL(GS, b)
response = b(1)./(b(2)+exp(b(3).*GS));
end