clear all
close all

%% Encoding Test
% See PAge 281 of Proakis Communication Systems Engineering
% Use this uncoded bitstream to create an encoded bitstream and a
% dictionary
print_flag = true; % Print flag for element by element comparison
uncoded_bitstream = ['0100001100001010000010100000110000010100001001001'];
N = size(uncoded_bitstream,2);
%uncoded_bitstream = reshape(uncoded_bitstream,[N,1]); % one long bitstream that is N x 1
num_bits = size(uncoded_bitstream,2);
bits_per_phrase = 4;   
[encoded_bitstream, encoding_scheme] = lempelzivencoding(uncoded_bitstream,num_bits,bits_per_phrase);
% Compression ratio of example
N_encoded = size(encoded_bitstream,2);
cr = N/N_encoded
actual_contents = ["","0","1","00","001","10","000","101","0000","01","010","00001","100","0001","0100","0010"];
% Print Dictionary versus actual from encoding scheme
for i = 0:15
    locations = encoding_scheme.locations;
    contents = encoding_scheme.contents;
    loc = i;
    con = contents(num2str(loc));
    if print_flag
        fprintf('loc: %u ==> content: %s  | book: %s \n', loc, con, actual_contents(i+1))
    end
    assert(con == actual_contents(i+1))
end
fprintf('Passed dictionary generation: Program dictionary is same as Pg 282 from textbook\n')

%% Decoding test
decoded_bitstream = lempelzivdecoding(encoded_bitstream,bits_per_phrase,encoding_scheme);
minlength = min(length(decoded_bitstream), length(uncoded_bitstream));
num_errors = 0;
for i = 1:minlength
    error  = (decoded_bitstream(i) ~= uncoded_bitstream(i));
    num_errors = num_errors + error;
    assert(num_errors == 0);
end

%% Further tests
% Use dictionary to encode the example from pg 281
% Stored as strings in order to hold them in an array 
% Convert strings into char in the test code
uncoded_ex = ["0";"1"; "00";"001";"10";"000";"101";"0000";"01";"010";"00001";"100";"0001";"0100";"0010";"01001"];
encoded_ex_actual = ["00000";"00001";"00010";"00111";"00100";"00110";"01011";"01100";"00011";"10010";"10001";"01010";"01101";"10100";"01000"; "11101"];
total_uncoded_length = 0;
total_encoded_length = 0;
for i = 1:16
    %ENCODED
    num_bits = strlength(uncoded_ex(i)); % number of bits inputted is strlength
    total_uncoded_length = total_uncoded_length + num_bits;
    [encoded_ex, ~]  = lempelzivencoding(convertStringsToChars(uncoded_ex(i)),num_bits,bits_per_phrase,encoding_scheme);
    total_encoded_length = total_encoded_length + strlength(encoded_ex);
    
    %decoded
    decoded_ex  = lempelzivdecoding(encoded_ex,bits_per_phrase,encoding_scheme);
    if print_flag
        fprintf('uncoded_ex: %s | encoded_ex: %s | encoded_ex_actual %s | decoded_ex: %s \n',uncoded_ex(i),encoded_ex,encoded_ex_actual(i), decoded_ex)
    end
    assert(encoded_ex_actual(i) == encoded_ex)
    assert(uncoded_ex(i) == decoded_ex);
end
cr_book = 49/80;
cr_ex = total_uncoded_length/total_encoded_length;
assert(abs(cr_book - cr_ex) < 1e-6) % asserts that the value are very close barring rounding precision
fprintf('Passed encoding confirmation: Program encoding and decoding matches pg 281 from textbook\n')

