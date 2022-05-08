% Function lempelziv
% Applies lempelziv encoding scheme
% Inputs:
%   - Uncoded bit stream as a character array of bits that is 1 x num_bits
%   - num_bits in bit stream
%   - bits_per_phase as the number of encoded bits per phrase ( without the
%   codeword
%   - encoding_scheme: If maps are already created, utilize the maps withing the encoding_scheme struct 
%Outputs:
%   - encoded_bitstream: char array of 1 x new_num_bits. Encoded using
%   lempel ziv with each phrase being bits_per_phrase + 1 bit codeword long
%   - encoding_scheme: as struct containing maps for location and contents.
function [encoded_bitstream, encoding_scheme] = lempelzivencoding (uncoded_bitstream,num_bits,bits_per_phrase, encoding_scheme)
    create_encoding = true; % If encoding_scheme is not defined, create encoding
    if exist('encoding_scheme','var')
        create_encoding = false;
        get_dict_locations = encoding_scheme.locations;
        get_dict_contents = encoding_scheme.contents;
    end
    if create_encoding
        % Develop dictionary from bit stream. key and input are both char array
        % Use Map, indexing with dictionary contents, to look up
            % dictionary location (i.e. find content 001 at location 3)
        % Another Map, indexing as dictionary location 0-(2^n-1)for n = num bits, to
            % look up dictionary location (i.e. location 3 has content 001)
            % dict_location of 0 corresponds with ''
        get_dict_locations = containers.Map; % Use content as field to get locaiton num
        get_dict_contents = containers.Map; % use location as field to get content
        % Fill DEfault values
        get_dict_locations('') = '0';
        get_dict_contents('0') = '';
        % Fill  Dictionary with 2^bits new sequences then quit
        content = '';
        location = 1;
        for bit = 1:num_bits
            content = strcat(content, uncoded_bitstream(bit)); % Append new bit to phrase
            if isKey(get_dict_locations,content) % if sequence is in dictionary already
                continue; % 
            elseif location < 2^bits_per_phrase % if sequence is not in dictionary already and its not full
                get_dict_locations(content) = int2str(location);
                get_dict_contents(int2str(location)) = content;
                location = location + 1;
                content = '';
            else % if sequence is not in dictionary and its full
                break;
            end
        end
    end
    
    %Store in encoding scheme
    encoding_scheme.locations = get_dict_locations;
    encoding_scheme.contents = get_dict_contents;

    % Encode the bit stream using the dictionary
    % (i.e. if 3-bit dictionary 3rd entry  is 001 and a unique
    % sequence is 0010 then the output is 011 0
    encoded_bitstream = [];
    content = '';
    bit = 1;
    while bit <= num_bits + 1
        if isKey(get_dict_locations,content) & (bit ~= num_bits+1) % If the content is in, append the current bit unless the last bit is appended already
            % Append current bit to content
            content = strcat(content,uncoded_bitstream(bit));
        else
           
            % Take the first part and then add the last bit as
            % codeword
            % Do one-bit scenario first
            if length(content) == 1 % if its '0' or '1'
                prev = ''; % prior
                codeword = content(end); % last 
            else
                prev = content(1:end-1);
                codeword = content(end);
            end


            %Encode by finding location of prev in dictionary and appending
            %codeword
            encoded_bits = '';
            %Quantizer is unsigned fixed with the respective number of bits
            % for indexing
            q = quantizer('ufixed', [bits_per_phrase 0]); 
            encoded_bits = strcat(encoded_bits,num2bin(q,str2num(get_dict_locations(prev))));
            encoded_bits = strcat(encoded_bits, codeword);
            encoded_bitstream = [encoded_bitstream,  encoded_bits];
            % Reset content 
            if bit <= num_bits
                content = uncoded_bitstream(bit);
            else
                content = '';
            end
            
        end
        bit = bit + 1;
    end
end