% Function audioencoding.m
% Encodes a series of samples with lempel ziv or other encoding schemes.
%Input:
%   input_sig_struct: struct containing input signal
%       - sig: array of doubles 
%       - Fs: sample frequency as a double
%   source_coding_type: string input classifying which source coding
%   bits_per_phrase: number of bits to encode into 
%Output:
%   encoded_bitstream: array of bits with audio encoding
%   encoding_scheme: struct containing
%       -get_dict_contents: dictionary to get content from a location
%       -get_dict_locations: dictionary to get locations from a content
%   coderate: number of uncoded bits/ encoded bits
function [encoded_bitstream, encoding_scheme, coderate, uncoded_bitstream] = audioencoding(input_sig_struct, source_coding_type, bits_per_phrase)
    sig = input_sig_struct.sig; % input signal  as N x 1 double array
    Fs = input_sig_struct.Fs;    
    
    switch source_coding_type
        case "LZ" %Lempel Ziv from pg 281 of CSE
            % Get Uncoded bit stream from signal
            q = quantizer('double');
            sig_bit = num2bin(q,sig); % input signal bit stream as N x 64 char array
            sig_bit = sig_bit.'; % 64 x N
            uncoded_bitsize = q.Format(1);
            N = size(sig,1);
            uncoded_bitstream = reshape(sig_bit,[1,N*uncoded_bitsize]); % one long bitstream that is 64*N x 1 taken sig_bit samples from the columns
            num_bits = length(uncoded_bitstream);
            
            [encoded_bitstream, encoding_scheme] = lempelzivencoding(uncoded_bitstream, num_bits, bits_per_phrase);
          
    end
    
    coderate = length(encoded_bitstream)/length(uncoded_bitstream);
end


