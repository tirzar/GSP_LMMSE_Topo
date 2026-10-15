function sigma2 = ut_snr(snrDB)
%UT_SNR  Noise variance from SNR [dB] = -10 log10(sigma^2) (sigma in p.u.).
sigma2 = 10 .^ (-snrDB / 10);
end
