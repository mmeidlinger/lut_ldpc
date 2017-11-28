function sig = snr2sig(snr, R)
    sig = 10.^(-snr/20)./sqrt(2*R);
end