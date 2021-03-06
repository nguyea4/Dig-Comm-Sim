% Function viterbidecoding utilizes hard viterbi decoding with certain parameters
% In program, yet to implement path memory truncation
% Inputs:
%   - received_bitstream: 1 x N char array of bits
%   - k: block length, number of bits being fed into shift register per stage (usually 1)
%   - L: constraint length, number of blocks in shift register
%   - n: n-linear combinations on shift register
%   - G: generator polynomial (assuming k=1)
% outputs:
%   - decoded_bitstream : 1 x M char array of bits
function decoded_bitstream = viterbidecoding(received_bitstream, k,L,n,G)
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

    next_state_table;
    output_table;

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

    input_table;

    % Initialize the path_metric
    % computationally viable approach
    q1 = quantizer('ufixed', [1,0]); % quantizer for integer 1 or 0 to '1' or '0'
    q = quantizer('ufixed', [L-1, 0]); % quantizer for 2 bits to get output
    path_metric_size = length(received_bitstream)/n;%min(5*L+1, length(received_bitstream)/n);  % if sequence is shorter than 5*L use the other one
    path_metric = -1*ones([2^(L-1), path_metric_size]); % [not implemented] PM stores 5L + 1 stages using path memory truncation
    seq_path = -1*ones([2^(L-1), path_metric_size]); % stores info on which node before it gave it hte smallest magnitude
    % Populate path_metric that shows the path metric for ending on that
    % specific state
    for i=0:n:length(received_bitstream)-2
        rec_num = received_bitstream(i+1:i+n); % int array of bits
        rec_word = num2bin(q1,rec_num); % char array of bits 1 x n
        rec_word = reshape(rec_word, [1,length(rec_word)]);
        if mod(i, 5*L+1) == 0
            % Start dynamic memory 
        end

        % i/n th stage, populate PM[s,i+1]
        for ns = 0:2^(L-1)-1 % iterate over next state
            ind0 = find(input_table(:,ns+1)==0); 
            ind1 = find(input_table(:,ns+1)==1); % Get index of current val
            temp  = [];
            temp_states = [];
            for in = 0:1
                ind = find(input_table(:,ns+1)==in); % Get index of current state that could yield the next state
                if size(ind,1) ~= 0  && size(ind,2) ~= 0
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
                        temp_states = [temp_states, ind(j)-1];
                    end
                end
            end
            % take the min
            [M,I] = min(temp);
            path_metric(ns+1,floor(i/n)+1) = M;
            seq_path(ns+1,floor(i/n)+1) = temp_states(I);
        end

    end
    path_metric;

    % Find the optimal path now
    %1) find state sequence
    state_sequence = [];
    decoded_bitstream = [];
    initial_state = 0;
    for i = path_metric_size:-1:1
        if i == path_metric_size % Pick whatever minimum, they are equiprobable
            [M,I] = min(path_metric(:,i)); % min M and row index I
        else
            I = seq_path(state_sequence(1)+1,i+1)+1; % to make it an index for previous sequence
            in = input_table(I, state_sequence(1)+1); 
            decoded_bitstream = [in, decoded_bitstream];
%             paths = find(input_table(:, state_sequence(1)+1)~=-1); % restrict to paths from previous sequence
%             [M,I] = min(path_metric(paths,i)); % min M and row index I of Current seq of those paths
%             %I = temp_states(state_sequence);
%             in = input_table(I, state_sequence(1)+1); % Get the input from (current seq I, last state sequence next seq)
%             decoded_bitstream = [in,decoded_bitstream];  % append input
        end
        state_sequence = [I-1, state_sequence]; % index ranges 0-3 not 1-4
        
    end
    in = input_table(initial_state+1, state_sequence(1)+1); % Put the first state in there
    decoded_bitstream = [in, decoded_bitstream];
    state_sequence = [initial_state, state_sequence];
    
    % Convert decoded bitstream bit int to bit char
    q = quantizer('ufixed', [1,0]);
    decoded_bitstream = num2bin(q,decoded_bitstream).'; % transpose to make it row vector
    
    % Remove k*(L-1) zero padding
    decoded_bitstream = decoded_bitstream(1:end-k*(L-1));
end