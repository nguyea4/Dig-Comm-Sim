function noisy_sig = awgn(input_sig, SNRdB)
    N = length(input_sig);
    SNRlin = 10^(SNRdB/10);
    sigpow = rms(input_sig)^2;
    noisepow = sigpow/SNRlin;
    noise = sqrt(noisepow)*randn(1,N);
    noisy_sig = input_sig+noise;
end