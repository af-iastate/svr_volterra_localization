clearvars; close all; clc;
% load Volterra filtering method scripts
currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
addpath(fullfile(pathstr, 'ofdm'));


% --- Init vars ---
fc = 25e6;
scs = 15e3;
nFFT = 128;
qam_bitsPerSymbol = 4;               % Bits per symbol
qam_nSymbols = 2^qam_bitsPerSymbol;  % QAM-16

tx_symbols = [0, 0, randi([0, qam_nSymbols-1], nFFT - 4, 1).', 0, 0].';
tx_grid = qammod(tx_symbols, qam_nSymbols, UnitAveragePower=true);

[tx_iq, Fs] = ofdm_modulate(tx_grid, scs);
x = iq_modulate(tx_iq, Fs, fc);


% --- Simulate Transmission Channel ---
y = awgn(x, 40); % Add noise


% --- Plot OFDM RS/TX signal ---
figure(1)
t = (0:length(x)-1)./Fs;
y_t = (0:length(y)-1)./Fs;

subplot(211)
plot(t, x)
hold on
plot(y_t, y)
legend({'TX', 'RX'})
title('RX/TX Time domain')

subplot(212)
pspectrum(x, Fs)
hold on
pspectrum(y, Fs)
legend({'TX', 'RX'})


% ---- Demodulate RX ---
headroom_coeff = 2;
rx_iq = iq_demodulate(y, Fs, fc, (headroom_coeff*scs*nFFT/2));
rx_grid = ofdm_demodulate(rx_iq, nFFT);


% --- Plot OFDM TX/RX Grid Compare ---
figure(2)
rx_t = (0:length(rx_iq)-1)./Fs;

subplot(211)
plot(t, real(tx_iq))
hold on
plot(rx_t, real(rx_iq))
legend({'TX', 'RX'})
title('Real Component')

subplot(212)
plot(t, imag(tx_iq))
hold on
plot(rx_t, imag(rx_iq))
legend({'TX', 'RX'})
title('Imaginary Component')


% --- Plot TX/RX Constellation ---
figure(3)
scatter(real(tx_grid), imag(tx_grid))
hold on
scatter(real(rx_grid), imag(rx_grid), 'x')
legend({'TX', 'RX'})
title('QAM-16 Constellation')


% --- Demodulate QAM-16 ---
rx_symbols = qamdemod(rx_grid, qam_nSymbols, UnitAveragePower=true);
if max(tx_symbols - rx_symbols) < 1e-8
    disp("Oversampled receiver output matches transmitter input.");
else
    disp("Received symbols do not match transmitted symbols.")
end
