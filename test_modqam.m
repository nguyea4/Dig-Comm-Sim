close all;
clear all;

% Test plot a few of them
% bitstream, Am, M, fc, symbol_period, fs
Am = 1; M = 16; fc = 10; symbol_period = 1; fs = 44000;
print_flag = false;
t = 0:1/fs:symbol_period;
test = ['0000';'0101';'0100';'0010'];
a = modqam(test(1,:), Am, M, fc, symbol_period, fs);
b = modqam(test(2,:),Am, M, fc, symbol_period, fs);
c = modqam(test(3,:),Am, M, fc, symbol_period, fs);
d = modqam(test(4,:),Am, M, fc, symbol_period, fs);
figure
plot(a); hold on; plot(b); plot(c); plot(d); hold off
title("Multiple codewords in parallel")
legend(test)

% Test multiple streams adjacent
e = modqam('100001000',Am,M,fc,symbol_period,fs);
figure
plot(e)
title("Codewords in series")

% Compare Acos(wt+ph) to Acos(wt)+Asin(wt)
a; % '0000' maps to [-1,1] on constellation
a_test = Am*(-1*cos(2*pi*fc*t)+1*sin(2*pi*fc*t));
figure
plot(a); hold on; plot(a_test); hold off;
legend('a','a IQ split')

% correlators, assume 0 carrier phase, fully synched
%a: -1, 1
%b = -1/3,1/3
psi1 = cos(2*pi*fc*t);
psi2 = sin(2*pi*fc*t);
const = 2/(Am*symbol_period); % Am/2 is rms?
corr_a1 = 1/fs*sum(const*psi1.*a);
corr_a2 = 1/fs*sum(const*psi2.*a);
corr_b1 = 1/fs*sum(const*psi1.*b);
corr_b2 = 1/fs*sum(const*psi2.*b);
fprintf("Expected -- A_I: %f, A_Q: %f, B_I: %f, B_Q: %f \n", -1, 1, -1/3, 1/3)
fprintf("Calculated -- A_I: %f, A_Q: %f, B_I: %f, B_Q: %f \n", corr_a1, corr_a2, corr_b1, corr_b2)
bitstream = demodqam(e, Am, M, fc, symbol_period, fs)
assert(strcmp(bitstream ,'100001000000')) % Adds 3 more 0 bits than in  e because we loaded it up to the right size

% Create SNR curves
if print_flag % Change this to run it
    Am = 1; M = 16; fc = 1000; symbol_period = .005; fs = 44000;
    SNRdBs = -15:1:-7;
    bit_num = 10^6;
    error_rate = ones(size(SNRdBs));
    % At different SNR
    for s = 1:length(SNRdBs)
        if s ~= 1 
            if error_rate(s-1) == 0
                break; % break if previous was 0 to stop from running
            end
        end
        % Create randomize bitstream of 10^-5
        bits = [];
        for i = 1:bit_num
            if round(rand) == 0
                bits = [bits,'0'];
            else
                bits = [bits,'1'];
            end
        end

        % modulated
        mod_sig = modqam(bits,  Am, M, fc, symbol_period, fs);
        
        %awgn
        noisy_sig = awgn(mod_sig, SNRdBs(s));
        
        %demodulate
        decoded_bits = demodqam(noisy_sig, Am, M, fc, symbol_period, fs);
        
        error = 0;
        for i = 1:length(decoded_bits)
            if decoded_bits(i) ~= bits(i)
                error = error + 1;
            end
        end
        
        error_rate(s) = error/length(decoded_bits);
        
    end
    figure
    semilogy(SNRdBs,error_rate);
    title("BER vs SNRdB")
end