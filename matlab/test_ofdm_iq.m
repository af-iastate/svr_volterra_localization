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

iq_fLo = 30e6;

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


% --- I/Q Modulate ---
% time vector
tx_tEnd = ofdm_TsUp * (length(tx_ofdmSymbolsUpFull)-1);
tx_t = (0:ofdm_TsUp:tx_tEnd).';
% i/q
tx_iqUp = real(tx_ofdmSymbolsUpFull).*cos(2*pi*iq_fLo.*tx_t) ...
            - imag(tx_ofdmSymbolsUpFull).*sin(2*pi*iq_fLo.*tx_t);

tx_output = tx_iqUp;

% --- Transmission Line Distort ---
rx_input = tx_output;
rx_input = awgn(rx_input, 40);  % Add noise

% --- Receiver I/Q Demodulate ---
rx_tEnd = ofdm_TsUp * (length(rx_input)-1);
rx_t = (0:ofdm_TsUp:rx_tEnd).';
[rx_iqFiltB, rx_iqFiltA] = butter(5, (2*ofdm_scs*ofdm_nFFT/2)/(ofdm_FsUp/2), 'low');
% i/q
rx_i = 2*filtfilt(rx_iqFiltB, rx_iqFiltA, rx_input .* cos(2*pi*iq_fLo.*rx_t));
rx_q = 2*filtfilt(rx_iqFiltB, rx_iqFiltA, rx_input .* -sin(2*pi*iq_fLo.*rx_t));
rx_iq = rx_i + 1j.*rx_q;

% --- Receiver Remove Offset ---
rx_sync = rx_iq(ofdm_nCycPreUp+1 : end);

% --- Plot Test ---
figure(2)
subplot(211)
plot(tx_t, real(tx_ofdmSymbolsUpFull))
hold on
plot(rx_t(ofdm_nCycPreUp+1 : end), real(rx_sync))
title('Real Component')
legend("TX", "RX", "Location", "southeast")

subplot(212)
hold on
plot(tx_t, imag(tx_ofdmSymbolsUpFull))
plot(rx_t(ofdm_nCycPreUp+1 : end), imag(rx_sync))
title('Imag Component')
legend("TX", "RX", "Location", "southeast")

% --- Receiver Un-OFDM ---
rx_gridUp = fft(rx_sync(1:ofdm_nFFTUp), ofdm_nFFTUp) / ofdm_upFactor;

% --- Receiver Equalize ---
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
