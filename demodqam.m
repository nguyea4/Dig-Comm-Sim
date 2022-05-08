function output_bitstream = demodqam(rec_sig, Am, M, fc, symbol_period, fs)
    % Assume synchronization
    % split into windows
    N = length(rec_sig);
    t = 0:1/fs:symbol_period; % for perfect window
    win_length = length(t);
    x = [-1*ones([1,4]), -1/3*ones([1,4]), 1*ones([1,4]), 1/3*ones([1,4])];
    temp = [1, 1/3, -1, -1/3];
    y = [temp, temp, temp, temp];
    constellation = x+y*1j;
    
    q = quantizer('ufixed', [4 0]); 
    bitstream = [];
    for i = 1:win_length:N
        if i+win_length -1 > N
            i  = N - win_length + 1; % to make it same legnth at all times
        end
        window = rec_sig(i:i+win_length-1);
        % Correlate with basis function between each fully synched window
        psi1 = cos(2*pi*fc*t);
        psi2 = sin(2*pi*fc*t);
        const_I = 1/fs*sum(2/(Am*symbol_period)*psi1.*window);
        const_Q = 1/fs*sum(2/(Am*symbol_period)*psi2.*window);
        constellation_est = const_I+1j*const_Q;
        
        % Compare I and Q to the respective modulationa nd choose shortest
        % euclidean distance. Find shortest distance
        min = -1; % temp holder
        nearest = 0;
        for j = 1:length(constellation)
            dist = abs(constellation_est - constellation(j));
            if min == -1 || dist < min
                nearest = constellation(j);
                min = dist;
            end
        end
        k = find(constellation == nearest)-1;
        bitstream = [bitstream,num2bin(q,k)];
    end
    output_bitstream = bitstream;
end