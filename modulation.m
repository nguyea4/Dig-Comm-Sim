%file hierarchy: modulation -> modqam -> constellation2signal
function [output_sig, fs]  = modulation(bitstream,fs, mod_scheme, parameters)
switch mod_scheme
    case 'QAM'
        M = parameters.M; % M-bit QAM
        fc = parameters.fc; % Carrier frequency
        Am = parameters.Am; % Amplitude
        symbol_period = parameters.symbol_period; % symbol period
        [output_sig, fs] = modqam(bitstream, Am, M, fc, symbol_period, fs);
        
end
end
