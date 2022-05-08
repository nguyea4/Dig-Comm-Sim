clear all
close all

%% Fill Parameters
filename = "FurElise.wav";
max_time = 1; % how long of an audio file do you want to use in seconds
bits_per_phrase = 12;

% Read Audio file  and create struct
[input_sig, Fs] = audioread(filename);
num_samples = size(input_sig, 1);
max_sample = round(Fs*max_time);% the maximum time sample to reach.
input_sig_struct.sig = input_sig(1:max_sample,1); %Keep only channel 1 for simplicity. 
input_sig_struct.Fs = Fs;

% Audioencoding - using lempel ziv
[encoded_bitstream,encoding_scheme, coderate, uncoded_bitstream] = audioencoding(input_sig_struct, "LZ", bits_per_phrase);
fprintf("Code rate: %f", coderate);

% Convolutional coding

% QAM modulation
mod_scheme = 'QAM';
parameters.M = 16; % M-bit QAM
parameters.fc = 3000; % Carrier frequency
parameters.Am = 1; % Amplitude
parameters.symbol_period = 500*1/Fs; % symbol period of 500 samples, at 44kHz, 10 ms symbol period
mod_sig = modulation(encoded_bitstream,Fs, mod_scheme, parameters);

% AWGN

% QAM demodulation


% Convolutional decoding

% 



