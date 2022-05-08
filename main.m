clear all
close all

%% Fill Parameters
filename = "FurElise_short.wav";
max_time = 1; % how long of an audio file do you want to use in seconds
bits_per_phrase = 12;
savefile = './tmp/SNR-10dB'

% Read Audio file  and create struct
[input_sig, Fs] = audioread(filename);
num_samples = size(input_sig, 1);
max_sample = round(Fs*max_time);% the maximum time sample to reach.
randstart = randi(num_samples-max_sample-10,1); % Create random sample in the music
input_sig_struct.sig = input_sig(randstart:randstart+max_sample,1); %Keep only channel 1 for simplicity. 
input_sig_struct.Fs = Fs;

% Audioencoding - using lempel ziv
[encoded_bitstream,encoding_scheme, compress_rate, uncoded_bitstream] = audioencoding(input_sig_struct, "LZ", bits_per_phrase);
fprintf("LZ Code rate: %f\n", compress_rate);
save(savefile)

% Convolutional coding
k = 1; % # of bits per block
L = 3; % Constraint length, number of blocks
n=2; % n-linear combination of shift register
g1 = [1 0 1]; % k x L
g2 = [1 1 1]; % k x L
G = [g1;g2]; % n x L when k = 1
[convcoded_bitstream, coderate] = convcode(encoded_bitstream,k,L,n,G);
fprintf("Conv code Rate: %f \n", coderate);
save(savefile)

% QAM modulation
mod_scheme = 'QAM';
parameters.M = 16; % M-bit QAM
parameters.fc = 3000; % Carrier frequency
parameters.Am = 1; % Amplitude
parameters.symbol_period = 500*1/Fs; % symbol period of 500 samples, at 44kHz, 10 ms symbol period
mod_sig = modulation(convcoded_bitstream,Fs, mod_scheme, parameters);
fprintf("QAM Modulation completed \n");
save(savefile)

% AWGN
SNRdB = -10;
noisy_sig = awgn(mod_sig,SNRdB); 
fprintf("Noise added \n");
save(savefile)

% QAM demodulation
decoded_bits = demodqam(noisy_sig, parameters.Am, parameters.M, parameters.fc, parameters.symbol_period, Fs);
fprintf("QAM demodulation completed \n");
save(savefile)

% Convolutional decoding
decoded_bitstream = viterbidecoding(decoded_bits, k,L,n,G);
fprintf("Viterbi decoding completed \n");
save(savefile)

% Uncompress received data
signal_bitstream = lempelzivdecoding(decoded_bitstream,bits_per_phrase, encoding_scheme);
fprintf(" Uncompressing lempel-ziv completed \n");
save(savefile)

% Bring back to signal waveform
q = quantizer('double');
num_bits = q.Format(1);
reconstructed_sig = [];
% make sure the signal is multiple of 64
signal_bitstream = signal_bitstream(1: floor(length(signal_bitstream)/64) * 64); % truncate length to be a multipl eof 64
for i = 1:num_bits:length(signal_bitstream)
    endindex = min(i+num_bits-1, length(signal_bitstream));
    signal_value = bin2num(q,signal_bitstream(i:endindex));
    reconstructed_sig = [reconstructed_sig, signal_value]; % input signal bit stream as N x 64 char array
end
audiowrite(strcat(savefile,'.wav'),reconstructed_sig,Fs);
save(savefile)


%% Calculate metrics
% QAM BER
min_length = min(length(convcoded_bitstream), length(decoded_bits));
qamerror = 0;
for i = 1: min_length
    qamerror = qamerror + (convcoded_bitstream(i) ~=decoded_bits(i));
end
qamerror = qamerror/min_length

% Convolutional coding and QAM BER, compare after LZ and after viterbi
% decoding
min_length = min(length(encoded_bitstream), length(decoded_bitstream));
converror = 0;
for i = 1: min_length
    converror = converror + (encoded_bitstream(i) ~=decoded_bitstream(i));
end
converror = converror/min_length
save(savefile)

% Reconstructed BER
min_length = min(length(uncoded_bitstream), length(signal_bitstream));
reconstructerror = 0;
for i = 1: min_length
    reconstructerror = reconstructerror + (uncoded_bitstream(i) ~=signal_bitstream(i));
end
reconstructerror = reconstructerror/min_length
save(savefile)

% MMSE
min_length = min(length(input_sig), length(reconstructed_sig));
mmse = 0;
for i = 1:min_length
    mmse = mmse + (input_sig(i)-reconstructed_sig(i)).^2;
end
mmse = 1/min_length * mmse




