clear all
close all

%% Fill Parameters

filename = "FurElise_short.wav";
max_time = 1; % how long of an audio file do you want to use in seconds
bits_per_phrase = 12;

% Read Audio file  and create struct
[input_sig, Fs] = audioread(filename);
num_samples = size(input_sig, 1);
max_sample = round(Fs*max_time);% the maximum time sample to reach.
randstart = randi(num_samples-max_sample-10,1); % Create random sample in the music
input_sig = input_sig(randstart:randstart+max_sample,1); %Keep only channel 1 for simplicity. 
input_sig_struct.sig = input_sig;
input_sig_struct.Fs = Fs;

% Audioencoding
[encoded_bitstream,encoding_scheme, coderate, uncoded_bitstream] = audioencoding(input_sig_struct, "LZ", bits_per_phrase);


% Audio decoding  from the encoded audio bit stream
[decoded_bitstream] = lempelzivdecoding(encoded_bitstream,bits_per_phrase, encoding_scheme);

%Compare decoded and uncoded
% Pick smallest size
decode_compare_size = length(decoded_bitstream);
if length(uncoded_bitstream) < length(decoded_bitstream)
    decode_compare_size = length(uncoded_bitstream);
end
num_errors = 0;
for i = 1:decode_compare_size
    bit_error = (uncoded_bitstream(i)~=decoded_bitstream(i));
    num_errors = num_errors + bit_error;
end
num_errors = num_errors 

% Bring back to signal waveform
signal_bitstream = decoded_bitstream;
q = quantizer('double');
num_bits = q.Format(1);
reconstructed_sig = [];
for i = 1:num_bits:length(signal_bitstream)
    endindex = min(i+num_bits-1, length(signal_bitstream));
    signal_value = bin2num(q,signal_bitstream(i:endindex));
    reconstructed_sig = [reconstructed_sig, signal_value]; % input signal bit stream as N x 64 char array
end

min_length = min(length(input_sig), length(reconstructed_sig));
mmse = 0;
for i = 1:min_length
    mmse = mmse + (input_sig(i)-reconstructed_sig(i)).^2;
    if mmse ~= 0
        continue;
    end
end
mmse = 1/min_length * mmse

audiowrite('./tmp/test.wav',reconstructed_sig, Fs)



