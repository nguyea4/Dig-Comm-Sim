% Assume perfect windowing
% bit stream as row vector
% output sig as row vector
function [output_sig, fs]  = modqam(bitstream, Am, M, fc, symbol_period, fs)
    % M is either 16 or 64
    % Normalize constellation to +1 real and + 1 imaginary
    if M == 16
        % 4 x 4 constellation: http://ecelabs.njit.edu/ece489v2/lab5.php
        % Split real and imag constellation to  -1 to 1 V 
        x = [-1*ones([1,4]), -1/3*ones([1,4]), 1*ones([1,4]), 1/3*ones([1,4])];
        temp = [1, 1/3, -1, -1/3];
        y = [temp, temp, temp, temp];
        points = x+y*j;
        
        % Get amplitude and phase
        A = abs(points);
        phase = angle(points);
        
        % Pad bitstream with zeros
        N = length(bitstream);
        if mod(N,4) ~=0
            new_N = N + (4-mod(N,4));
            for i = 1:new_N-N
                bitstream = [bitstream, '0'];
            end
            N = new_N;
            
        end
        
        %Modulation
        output_sig = [];
        t = 0:1/fs:symbol_period; % for perfect window
        q = quantizer('ufixed', [4 0]); 
        % Begin modulating 16 at a time
        for i = 1:4:N % i is start of each 16-bits
            phrase_num = bin2num(q,bitstream(i:i+3));
            A_i = A(phrase_num+1);
            phase_i = phase(phrase_num+1);
            window = Am*A_i*cos(2*pi*fc*t-phase_i); % Equivalent to inphase and quadrature portion
            output_sig = [output_sig, window];
        end
        
    elseif M == 64
        assert(false, "Not implemented");
    end
    
    % Get constellation
    

end