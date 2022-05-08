%% This file is used to test and understand convolutional coding
clear all;
close all;

input = ['1101011'];
exp_output = ['111010000100101011'];

% convert to numerical array
input_num = input - '0';
exp_output_num = exp_output - '0';

% Define properties
k = 1; % # of bits per block
L = 3; % Constraint length, number of blocks
n=2; % n-linear combination of shift register
Rc = k/n;
g1 = [1 0 1]; % k x L
g2 = [1 1 1]; % k x L
G = [g1;g2]; % n x L when k = 1

%% Encoding 
% After receiving bit sequence, run through trellis and find the most
% likely path
output  = [];

% Register as a stack of k x L
shiftreg = zeros([L,k]); % bottom is FI and top is FO

% Append k(L-1) zeros
input_num = [input_num, zeros([1,k*(L-1)])];

% Start pushing the shift stuff until last appended 0 is in
for i = 1:k:length(input_num)
    shiftreg(1:end-1,:) = shiftreg(2:end,:); % Move new shift register up
    shiftreg(end,:) = input_num(i:i+k-1); % feed in next k bits to SR
    
    % compute n new bits
    n1 = mod(sum(G(1,:)*shiftreg),2); % sum (k x L * k x L)  to 1 digit
    n2 = mod(sum(G(2,:)*shiftreg),2);
    output = [output, n1, n2];
end

output = (output)
exp_output_num = (exp_output_num)
assert(isequal(output, exp_output_num))

[output_bitstream, coderate] = convcode(input,k,L,n,G); % Testing homemade functions
assert(isequal(output_bitstream, exp_output));
fprintf("Convolutional code for k = %d, L = %d, from pg 624 worked successfully\n", k, L);

%% Decoding 
% Dyanmic programming: http://web.mit.edu/6.02/www/f2011/handouts/8.pdf
% Data structure: https://faraday.emu.edu.tr/eaince/ee562/Public/tutorial_ce_va2.htm
k = 1; % # of bits per block
L = L; % Constraint length, number of blocks
n=n; % n-linear combination of shift register
received_bitstream = output;

% Create various data structures

% have a next state table and output table
% Rows correspond to L-1 bit word like '00' for 2 bit, columns are for inputs of k-bits
% for L = 3, each column is '00','01,.... and k=1, columns are input '0'
% and '1'
% Ex: Shifting '1' into '01' -> '10'. index - 1is value of ph4
next_state_table = zeros([2^(L-1), 2^k]); 
output_table = zeros([2^(L-1),2^k]);
q = quantizer('ufixed', [L, 0]); % Quantizer to get 3 bits for generating poly
for cs = 0:2^(L-1)-1
    for in = 0:1
        % For next state table
        % bit shift current state(cs) to the right one '10' -> '01'
        % then add new bit, in, by adding bit*2^(L-2), if bit = 1 add 2 to
        % get '11' in this example
        cs_ind = cs +1;
        in_ind = in + 1;
        next_state_table(cs_ind, in_ind) = bitshift(cs,-1) + in*2^(L-2);
        
        % For output table, utilize the generating polynomial
        out = 0;
        bit_arr = num2bin(q,cs+in*2^(L-1))-'0'; %0 -> '000' -> 0 0 0
        for j = 1:n
            nj = mod(sum(G(j,:)*bit_arr.'),2); % sum (k x L * k x L)  to 1 digit
            out = out + nj*2^(n-j);
        end
        output_table(cs_ind, in_ind) = out;
    end
end

next_state_table
output_table

input_table = -1*ones([2^(L-1), 2^(L-1)]); % For the current state(row) and next state (col) what input can obtain this

for cs = 0:2^(L-1)-1
    for ns = 0:2^(L-1)-1
        cs_ind = cs+1;
        ns_ind = ns + 1; % for matlab indexing
        ind0 = find(next_state_table(cs_ind,:) == ns);
        if ind0 % Has a value
            input_table(cs_ind, ns_ind) = ind0 -1; % put respective single bit(assume k = 1)
        end
    end
end

input_table

% Initialize the path_metric and optimal_prev
% computationally viable approach
q1 = quantizer('ufixed', [1,0]); % quantizer for integer 1 or 0 to '1' or '0'
q = quantizer('ufixed', [L-1, 0]); % quantizer for 2 bits to get output
initial_state = 0;
path_metric_size = min(5*L+1, length(received_bitstream)/n);  % if sequence is shorter than 5*L use the other one
path_metric = -1*ones([2^(L-1), path_metric_size]); % PM stores 5L + 1 stages using path memory truncation
optimal_prev = -1*ones([2^(L-1), path_metric_size]); %OP adjusts the optimal prev at each step to traceback
% Populate path_metric that shows the path metric for ending on that
% specific state
for i=0:n:length(received_bitstream)-2
    rec_num = received_bitstream(i+1:i+n); % int array of bits
    rec_word = num2bin(q1,rec_num).'; % char array of bits 1 x n
    if mod(i, 5*L+1) == 0
        % Start dynamic memory 
    end
    
    % i/n th stage, populate PM[s,i+1]
    for ns = 0:2^(L-1)-1 % iterate over next state
        ind0 = find(input_table(:,ns+1)==0); 
        ind1 = find(input_table(:,ns+1)==1); % Get index of current val
        temp  = [];
        for in = 0:1
            ind = find(input_table(:,ns+1)==in); % Get index of current state that could yield the next state
            for j = 1:length(ind) % For bits 0
                % if first stage, past hamming is 0, ow populate it
                if i == 0
                    past = 0;
                else
                    past = path_metric(ind(j),floor((i-1)/n)+1);
                end

                out = output_table(ind(j),in+1); % Output state as int
                out_bit = num2bin(q,out); % output as char array
                t  = abs(out_bit - rec_word); % '0' - '1' = -1 so abs to get 1
                t1 = sum(t); % number of bit flipped between the two
                t2 = t1 + past;
                temp = [temp, t2];
            end
        end
        % take the min
        path_metric(ns+1,floor(i/n)+1) = min(temp);
    end
    
end
path_metric

% Find the optimal path now
%1) find state sequence
state_sequence = [];
decoded_input = [];
initial_state = 0;
for i = path_metric_size:-1:1
    if i == path_metric_size % Pick whatever minimum, they are equiprobable
        [M,I] = min(path_metric(:,i)); % min M and row index I
    else
        paths = find(input_table(:, state_sequence(1)+1)~=-1); % restrict to paths from previous sequence
        [M,I] = min(path_metric(paths,i)); % min M and row index I of Current seq of those paths
        I = paths(I);
        in = input_table(I, state_sequence(1)+1); % Get the input from (current seq I, last state sequence next seq)
        decoded_input = [in,decoded_input];  % append input
    end
    state_sequence = [I-1, state_sequence]; % index ranges 0-3 not 1-4
end
in = input_table(initial_state+1, state_sequence(1)+1); % Put the first state in there
decoded_input = [in, decoded_input]
state_sequence = [initial_state, state_sequence]
decoded_bitstream = viterbidecoding(received_bitstream, k,L,n,G);
assert(isequal(decoded_input, input_num));
assert(isequal(decoded_bitstream,input));
fprintf("Viterbi Hard decoding check complete from example in textbook\n");


