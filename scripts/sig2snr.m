function snr = sig2snr(sig, R)
    snr = -10*log10(2*R.*sig.^2);
end