% Function audiodecoding.m
% Decodes a series of samples with lempel ziv or other encoding schemes.
%Input:
%   received_sig_struct: struct containing received signal
%       - sig: array of doubles 
%       - Fs: sample frequency as a double
%   encoding_scheme: struct containing mappin
%       -get_dict_contents: dictionary to get content from a location
%       -get_dict_locations: dictionary to get locations from a content
%   bits_per_phrase: number of bits to decode from (excluding codeword)
%Output:
%   decoded_bitstream: array of bits with audio encoding
%   decoding_scheme: struct containing
%       -get_dict_contents: dictionary to get content from a location
%       -get_dict_locations: dictionary to get locations from a content
%   coderate: number of uncoded bits/ encoded bits
function [decoded_bitstream] = audiodecoding(received_sig_struct, encoding_scheme, bits_per_phrase)
    sig = received_sig_struct.sig; % input signal  as N x 1 double array
    Fs = received_sig_struct.Fs;    
    
    switch source_coding_type
        case "LZ" %Lempel Ziv from pg 281 of CSE
            % Get encoded bit stream from received signal
            q = quantizer('double');
            sig_bit = num2bin(q,sig); % input signal bit stream as N x 64 char array
            encoded_bitsize = q.Format(1);
            N = size(sig,1);
            encoded_bistream = reshape(sig_bit,[N*encoded_bitsize,1]); % one long bitstream that is 64*N x 1
            num_bits = size(encoded_bistream,1);
            
            % Convert
            [decoded_bitstream] = lempelzivdecoding(encoded_bistream, bits_per_phrase, encoding_scheme);
    end
end

