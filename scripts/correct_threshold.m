%% Correct the SNR (E_b / N0) if it has been calculated based on a wrong rate.
% Input
%   snr1 ... Wrong SNR
%   R1   ... Wrong Rate
%   R2   ... True Rate
% Output
%   snr2 ... True SNR

function snr2 =  correct_snr(snr1, R1, R2)
    snr2 = snr1 + 10*log10(R1/R2);
end