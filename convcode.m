% Function convcode utilizes convolutional encoding with certain parameters
% Inputs:
%   - input_bitstream: 1 x N char array of bits
%   - k: block length, number of bits being fed into shift register per stage (usually 1)
%   - L: constraint length, number of blocks in shift register
%   - n: n-linear combinations on shift register
%   - G: generator polynomial (assuming k=1)
% outputs:
%   - output_bitstream: 1 x M char array of bits
%   - coderate: k/n, number of message bits per total bits
function [output_bitstream, coderate] = convcode(input_bitstream,k,L,n,G)
output_bitstream  = [];

coderate = k/n;

% Register as a stack of k x L
shiftreg = zeros([L,k]); % bottom is FI and top is FO

% Append k(L-1) zeros
input_bitstream = [input_bitstream, zeros([1,k*(L-1)])];

% Start pushing the shift stuff until last appended 0 is in
for i = 1:k:length(input_bitstream)
    shiftreg(1:end-1,:) = shiftreg(2:end,:); % Move new shift register up
    shiftreg(end,:) = input_bitstream(i:i+k-1); % feed in next k bits to SR
    
    % compute n new bits
    n1 = mod(sum(G(1,:)*shiftreg),2); % sum (k x L * k x L)  to 1 digit
    n2 = mod(sum(G(2,:)*shiftreg),2);
    output_bitstream = [output_bitstream, n1, n2];
end

% Convert ouput_bitstream as int array to char array
q = quantizer('ufixed', [1,0]);
output_bitstream = num2bin(q,output_bitstream).';

end