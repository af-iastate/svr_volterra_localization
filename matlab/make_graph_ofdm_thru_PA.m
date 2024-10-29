clearvars; close all; clc;
% load Volterra filtering method scripts
currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
addpath(fullfile(pathstr, 'ofdm'));
addpath(fullfile(pathstr, 'known-inverse'));
addpath(fullfile(pathstr, 'methods/morhac'));
addpath(fullfile(pathstr, 'methods/pm'));
addpath(fullfile(pathstr, 'methods/ppm'));


% --- Init OFDM ---
fc = 60e6;
scs = 15e3;
nFFT = 128;
qam_bitsPerSymbol = 4;               % Bits per symbol
qam_nSymbols = 2^qam_bitsPerSymbol;  % QAM-16

tx_symbols = [0, 0, randi([0, qam_nSymbols-1], nFFT - 4, 1).', 0, 0].';
tx_grid = qammod(tx_symbols, qam_nSymbols, UnitAveragePower=true);

[tx_iq, Fs] = ofdm_modulate(tx_grid, scs);
x = iq_modulate(tx_iq, Fs, fc);


% --- Init Predistorter ---
% system def
memblock_a = 0.999*GenBandstopIr(17, Fs, 25e6, 40e6);
polblock_b = 0.9*[1, 0.0, 0.2];
memblock_c = [0.98, 0.1, -0.3, 0.2];

% Build known-inverse kernel
M = length(memblock_a); 
K = 1; 
[h_true, K] = ApplyPolymap(memblock_a, M, K, polblock_b);

% Build forward system
polblock_b_f = @(x) ThirdOrderActivationFunc(x, polblock_b(1), polblock_b(3));
memblock_a_iir_b = [1, zeros(1, length(memblock_a) - 1)];
memblock_a_iir_a = memblock_a(:)';
memblock_a_f = @(x) filter(memblock_a_iir_b, memblock_a_iir_a, x);

forward_f = @(x) 0.88*memblock_a_f(polblock_b_f(x));


% --- Simulate Transmission Channel ---
y = forward_f(x);


% ---- Demodulate RX ---
headroom_coeff = 1.2;
rx_iq = iq_demodulate(y, Fs, fc, (headroom_coeff*scs*nFFT/2));
rx_grid = ofdm_demodulate(rx_iq, nFFT);


% --- Plot TX/RX Constellation ---
figure(3)
scatter(real(tx_grid), imag(tx_grid))
hold on
scatter(real(rx_grid), imag(rx_grid), 'x')
l_h = legend({'$TX$', '$RX$'}, ...
    'location', 'north', ...
    'Orientation', 'vertical', ...
    'Interpreter', 'latex');
l_h.Position = [0.659314697023859,0.723369380341809,0.228283777776378,0.103778507173644];
title('QAM-16 Constellation', 'Interpreter', 'latex')
xlabel('$I$', 'Interpreter', 'latex')
ylabel('$Q$', 'Interpreter', 'latex')
pbaspect([1,1,1])
drawnow

% --- Demodulate QAM-16 ---
rx_symbols = qamdemod(rx_grid, qam_nSymbols, UnitAveragePower=true);
if max(tx_symbols - rx_symbols) < 1e-8
    disp("Oversampled receiver output matches transmitter input.");
else
    disp("Received symbols do not match transmitted symbols.")
end


EVMs = abs(tx_grid - rx_grid) ./ abs(tx_grid);

avg_EVM = mean(EVMs);
worst_EVM = max(EVMs);

m2db = @(x) 20.*log10(x);

avg_EVM_dB = m2db(avg_EVM);
worst_EVM_dB = m2db(worst_EVM);

disp(avg_EVM_dB)
disp(worst_EVM_dB)
