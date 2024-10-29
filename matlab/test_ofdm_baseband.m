clc; close all; clearvars;

% --- Init Vars ---
qam_bitsPerSymbol = 4;               % Bits per symbol
qam_nSymbols = 2^qam_bitsPerSymbol;  % QAM-16

ofdm_nFFT = 128;                    % Number of OFDM FFT bins
ofdm_nCycPre = 32;                  % Number of Cyclic Prepend Indices
ofdm_scs = 15e3;                    % Subcarrier spacing in Hz
ofdm_Fs = ofdm_scs * ofdm_nFFT/2;   % Base Sampling rate 
ofdm_Ts = 1 / ofdm_Fs;              % Base Sample duration in seconds
ofdm_upFactor = 131;                % Up-sample factor
ofdm_FsUp = ofdm_upFactor* ofdm_Fs;
ofdm_TsUp = 1 / ofdm_FsUp;

% --- Generate Input Symbols ---
tx_symbols = randi([0, qam_nSymbols-1], ofdm_nFFT, 1);

% --- Generate QAM Grid ---
% base
tx_grid = qammod(tx_symbols, qam_nSymbols, UnitAveragePower=true);
% oversampled
tx_gridUp = [
    tx_grid(1 : ofdm_nFFT/2);
    zeros((ofdm_upFactor-1) * ofdm_nFFT, 1);
    tx_grid((ofdm_nFFT/2 +1) : ofdm_nFFT)
];

% --- Generate OFDM time-domain ---
% base
tx_ofdmSymbols = ifft(tx_grid, ofdm_nFFT);
tx_ofdmSymbols = tx_ofdmSymbols(:);
% oversampled
ofdm_nFFTUp = ofdm_upFactor * ofdm_nFFT;
tx_ofdmSymbolsUp = ofdm_upFactor * ifft(tx_gridUp, ofdm_nFFTUp);
tx_ofdmSymbolsUp = tx_ofdmSymbolsUp(:);

% --- Prepend Cyclical Prefix ---
% base
tx_ofdmSymbolsPrefix = tx_ofdmSymbols(ofdm_nFFT-ofdm_nCycPre+1 : ofdm_nFFT);
tx_ofdmSymbolsFull = [tx_ofdmSymbolsPrefix; tx_ofdmSymbols];
% oversampled
ofdm_nCycPreUp = ofdm_nCycPre * ofdm_upFactor;
tx_ofdmSymbolsUpPrefix = tx_ofdmSymbolsUp(ofdm_nFFTUp-ofdm_nCycPreUp+1 : ofdm_nFFTUp);
tx_ofdmSymbolsUpFull = [tx_ofdmSymbolsUpPrefix; tx_ofdmSymbolsUp];
tx_output = tx_ofdmSymbolsUpFull;

% --- Plot ---
figure(1)
% base
Tend = ofdm_Ts * (length(tx_ofdmSymbolsFull)-1);
subplot(211)
hold off
plot(0:ofdm_Ts:Tend, real(tx_ofdmSymbolsFull), "*")
title("Real component of transmitter output")
subplot(212)
hold off
plot(0:ofdm_Ts:Tend, imag(tx_ofdmSymbolsFull), "*")
title("Imaginary component of transmitter output")
% oversampled
TendUp = ofdm_TsUp * (length(tx_ofdmSymbolsUpFull)-1);
subplot(211)
hold on
plot(0:ofdm_TsUp:TendUp, real(tx_ofdmSymbolsUpFull))
legend("Original", "Upsampled", "Location", "southeast")
subplot(212)
hold on
plot(0:ofdm_TsUp:TendUp, imag(tx_ofdmSymbolsUpFull))
legend("Original", "Upsampled", "Location", "southeast")

% --- Transmission Line Dist ---
h_chan = [0.4 1 0.4].';
h_offset = (randi(ofdm_nCycPreUp) - 1);
channelDelay = dsp.Delay(1);

rx_input = tx_output;
rx_input = awgn(rx_input, 40);  % Add noise
% rx_input = conv(rx_input, h_chan);  % Apply Channel Filter
% rx_input = channelDelay(rx_input);  % Add Channel Delay
% rx_input = rx_input(ofdm_nCycPreUp+1+channelDelay.Length-h_offset : end);  % Add Channel Offset

% --- Receiver Remove Offset ---
rx_sync = rx_input(ofdm_nCycPreUp+1 : end);
rx_gridUp = fft(rx_sync(1:ofdm_nFFTUp), ofdm_nFFTUp) / ofdm_upFactor;

% --- Receiver Equalize ---
% H_chan = fft(h_chan, ofdm_nFFTUp);
% % Linear phase term related to timing offset
% H_offset = exp(-1j * 2*pi*h_offset * (0:ofdm_nFFTUp-1).'/ofdm_nFFTUp);
% rx_gridUpEq = rx_gridUp ./ (H_chan .* H_offset);
rx_gridUpEq = rx_gridUp;

% --- Receiver Downsample ---
rx_gridEq = [
    rx_gridUpEq(1: ofdm_nFFT/2);
    rx_gridUpEq((1+(ofdm_upFactor-0.5)*ofdm_nFFT) : end)
];

% --- Plot Recieved constellation ---
figure(3)
scatter(real(tx_grid), imag(tx_grid))
hold on
scatter(real(rx_gridEq), imag(rx_gridEq), 'x')

% --- Demodulate QAM-16 ---
rx_symbols = qamdemod(rx_gridEq, qam_nSymbols, UnitAveragePower=true);
if max(tx_symbols - rx_symbols) < 1e-8
    disp("Oversampled receiver output matches transmitter input.");
else
    disp("Received symbols do not match transmitted symbols.")
end
