% Function lempelziv
% Applies lempelziv decoded scheme
% Inputs:
%   - encoded bit stream as a character array of bits that is 1 x num_bits
%   - bits_per_phase as the number of encoded bits per phrase ( without the
%   codeword
%   - encoding_scheme:  utilize the maps withing the encoding_scheme struct 
%Outputs:
%   - decoded_bitstream: char array of 1 x new_num_bits. DEcode  using
%   lempel ziv with each phrase being bits_per_phrase + 1 bit codeword long
function [decoded_bitstream] = lempelzivdecoding (encoded_bit_stream,bits_per_phrase, encoding_scheme)
    get_dict_locations = encoding_scheme.locations;
    get_dict_contents = encoding_scheme.contents;
    
    num_bits = length(encoded_bit_stream);

    % Decode the bit stream using the dictionary
    % (i.e. if bits_per_phase = 4, we take sequence of 5 bits like 0010 0
    % and search the 0010=2nd location for the content at that and add a 0)
    decoded_bitstream = [];
    phrase = '';
    for bit = 1:bits_per_phrase+1:num_bits % should always have the same amount
        %Split phrase into location + codebit
        phrase = encoded_bit_stream(bit:bit+bits_per_phrase); % Get the bits_per_phrase+1 bits
        loc = phrase(1:end-1);
        codeword = phrase(end);
        
        % convert location from bit char array like '0000' to char of integer like '0'
        q = quantizer('ufixed', [bits_per_phrase 0]); 
        loc_num = bin2num(q,loc);
        loc = int2str(loc_num);
        
        % Get contents at that location and add that + codeword
        decoded_bits = '';
        decoded_bits = strcat(decoded_bits, get_dict_contents(loc),codeword);
        decoded_bitstream = [decoded_bitstream, decoded_bits];
    end
end