clearvars; close all; clc;
% load Volterra filtering method scripts
currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
addpath(fullfile(pathstr, 'ofdm'));
addpath(fullfile(pathstr, 'known-inverse'));
addpath(fullfile(pathstr, 'methods/direct'));
addpath(fullfile(pathstr, 'methods/morhac'));
addpath(fullfile(pathstr, 'methods/pm'));
addpath(fullfile(pathstr, 'methods/ppm'));
addpath(fullfile(pathstr, 'predist'));


% --- Init OFDM ---
fc = 20e6;
scs = 15e3;
nFFT = 2^7;
nCycPre = 32;
upfactor = 120;
qam_bitsPerSymbol = 4;               % Bits per symbol
qam_nSymbols = 2^qam_bitsPerSymbol;  % QAM-16
headroom_coeff = 1.5;


% --- Init Volterra RLS ---
rls_factor = upfactor / 8;
K_predist = 3;
M_predist = 7; %* (upfactor / rls_factor);
lambda = 1;

redo_sim = false;
sim_filename = 'volterra_iq_train.mat';

f_d = @(x) downsample(x, rls_factor);
f_u = @(x) sinc_upsample(x, rls_factor);


% --- Init Training ---
rng(20)
tx_train_symbols = [0, 0, randi([0, qam_nSymbols-1], nFFT - 4, 1).', 0, 0].';
tx_train_grid = qammod(tx_train_symbols, qam_nSymbols, UnitAveragePower=true);

[x_train_iq, Fs] = ofdm_modulate(tx_train_grid, scs, 'nCycPre', nCycPre, 'upfactor', upfactor);
x_train_iq_ds = f_d(x_train_iq);
x_train = iq_modulate(x_train_iq, Fs, fc);


% --- Init Testing ---
rng(29)
tx_test_symbols = [0, 0, randi([0, qam_nSymbols-1], nFFT - 4, 1).', 0, 0].';
tx_test_grid = qammod(tx_test_symbols, qam_nSymbols, UnitAveragePower=true);

[x_test_iq, ~] = ofdm_modulate(tx_test_grid, scs, 'nCycPre', nCycPre, 'upfactor', upfactor);
x_test_iq_ds = f_d(x_test_iq);
x_test = iq_modulate(x_test_iq, Fs, fc);


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

forward_f = @(x) memblock_a_f(polblock_b_f(x));


% --- Simulate Transmission Channel ---
y_train = forward_f(x_train);
y_train_iq = iq_demodulate(y_train, Fs, fc, (headroom_coeff*scs*nFFT/2));
y_train_iq_ds = f_d(y_train_iq);


% --- Plot Training OFDM RS/TX signal ---
figure(1)
t = (0:length(x_train)-1)./Fs;
y_t = (0:length(y_train)-1)./Fs;

subplot(211)
plot(t, x_train)
hold on
plot(y_t, y_train)
legend({'x', 'y'})
title('Training (No predistortion)')

subplot(212)
pspectrum(x_train, Fs)
hold on
pspectrum(y_train, Fs)
legend({'x', 'y'})


% --- Plot Training IQ RS/TX signal ---
figure(2)
t = (0:length(x_train_iq)-1)./Fs;
y_t = (0:length(y_train_iq)-1)./Fs;

subplot(211)
plot(t, real(x_train_iq))
hold on
plot(y_t, real(y_train_iq))
legend({'x(t)', 'y(t)'})
title('Training I/Q (No Predistortion) Real Component')

subplot(212)
plot(t, imag(x_train_iq))
hold on
plot(y_t, imag(y_train_iq))
legend({'x(t)', 'y(t)'})
title('Training I/Q (No Predistortion) Imaginary Component')


% --- Train SVR models ---
offset = 1 + nCycPre*(upfactor / rls_factor);
nSamples = nFFT * (upfactor / rls_factor);

if (exist(sim_filename, 'file') && ~redo_sim)
    load(sim_filename, 'h_i', 'P_i', 'h_q', 'P_q')
else
    [h_i, P_i, h_q, P_q] = volterraRlsOfflineIq(y_train_iq_ds, ...
                                            x_train_iq_ds, ...
                                            M_predist, ...
                                            K_predist, ...
                                            'lambda', lambda);
    save(sim_filename, 'h_i', 'P_i', 'h_q', 'P_q');
end


% --- ReSimulate Training Data with Predist ---
u_train_iq = f_u(volterraFilterIqDirect(x_train_iq_ds, h_i, M_predist, K_predist)) ...
    + 1j*f_u(volterraFilterIqDirect(x_train_iq_ds, h_q, M_predist, K_predist));
u_train = iq_modulate(u_train_iq, Fs, fc);
y_train_hat = forward_f(u_train);
y_train_hat_iq = iq_demodulate(y_train_hat, Fs, fc, (headroom_coeff*scs*nFFT/2));


% --- Plot Training OFDM RS/TX with Predistortion ---
figure(3)
t = (0:length(u_train)-1)./Fs;
y_t = (0:length(y_train_hat)-1)./Fs;

subplot(211)
plot(t, x_train)
hold on
plot(t, u_train)
plot(y_t, y_train_hat)
legend({'$x$', '$u$', '$\hat{y}$'}, 'Interpreter', 'latex')
title('Training (Predistortion)')

subplot(212)
pspectrum(x_train, Fs)
hold on
pspectrum(u_train, Fs)
pspectrum(y_train_hat, Fs)
legend({'$x$', '$u$', '$\hat{y}$'}, 'Interpreter', 'latex')


% --- Plot Training IQ RS/TX signal with Predistortion ---
figure(4)
t = (0:length(x_train_iq)-1)./Fs;
y_t = (0:length(y_train_hat_iq)-1)./Fs;

subplot(211)
plot(t, real(x_train_iq))
hold on
plot(y_t, real(y_train_hat_iq))
legend({'$x$',  '$\hat{y}$'}, 'Interpreter', 'latex')
title('Training I/Q (Predistortion) Real Component')

subplot(212)
plot(t, imag(x_train_iq))
hold on
plot(y_t, imag(y_train_hat_iq))
legend({'$x$', '$\hat{y}$'}, 'Interpreter', 'latex')
title('Training I/Q (Predistortion) Imaginary Component')


% ---- Demodulate Training RX ---
rx_train_iq = y_train_hat_iq;
rx_train_grid = ofdm_demodulate(rx_train_iq, nFFT, 'nCycPre', nCycPre, 'upfactor', upfactor);


% --- Plot Training TX/RX Constellation ---
figure(6)
scatter(real(tx_train_grid), imag(tx_train_grid))
hold on
scatter(real(rx_train_grid), imag(rx_train_grid), 'x')
legend({'TX', 'RX'})
title('Training QAM-16 Constellation')


% --- Demodulate Training QAM-16 ---
rx_train_symbols = qamdemod(rx_train_grid, qam_nSymbols, UnitAveragePower=true);
if max(tx_train_symbols - rx_train_symbols) < 1e-8
    disp("Oversampled Training receiver output matches transmitter input.");
else
    disp("Received Training symbols do not match transmitted symbols.")
end


% --- Simulate Test Data with Predist ---
u_test_iq = f_u(volterraFilterIqDirect(x_test_iq_ds, h_i, M_predist, K_predist)) ...
    + 1j*f_u(volterraFilterIqDirect(x_test_iq_ds, h_q, M_predist, K_predist));
u_test = iq_modulate(u_test_iq, Fs, fc);
y_test_hat = forward_f(u_test);
y_test_hat_iq = iq_demodulate(y_test_hat, Fs, fc, (headroom_coeff*scs*nFFT/2));


% --- Plot Training OFDM RS/TX with Predistortion ---
figure(7)
t = (0:length(u_test)-1)./Fs;
y_t = (0:length(y_test_hat)-1)./Fs;

subplot(211)
plot(t, x_test)
hold on
plot(t, u_test)
plot(y_t, y_test_hat)
legend({'$x$', '$u$', '$\hat{y}$'}, 'Interpreter', 'latex')
title('Test (Predistortion)')

subplot(212)
pspectrum(x_test, Fs)
hold on
pspectrum(u_test, Fs)
pspectrum(y_test_hat, Fs)
legend({'$x$', '$u$', '$\hat{y}$'}, 'Interpreter', 'latex')


% --- Plot Training IQ RS/TX signal with Predistortion ---
figure(8)
t = (0:length(x_test_iq)-1)./Fs;
y_t = (0:length(y_test_hat_iq)-1)./Fs;

subplot(211)
plot(t, real(x_test_iq))
hold on
plot(y_t, real(y_test_hat_iq))
legend({'$x$',  '$\hat{y}$'}, 'Interpreter', 'latex')
title('Test I/Q (Predistortion) Real Component')

subplot(212)
plot(t, imag(x_test_iq))
hold on
plot(y_t, imag(y_test_hat_iq))
legend({'$x$', '$\hat{y}$'}, 'Interpreter', 'latex')
title('Test I/Q (Predistortion) Imaginary Component')


% ---- Demodulate Training RX ---
rx_test_iq = y_test_hat_iq;
rx_test_grid = ofdm_demodulate(rx_test_iq, nFFT, 'nCycPre', nCycPre, 'upfactor', upfactor);


% --- Plot Training TX/RX Constellation ---
figure(9)
scatter(real(tx_test_grid), imag(tx_test_grid))
hold on
scatter(real(rx_test_grid), imag(rx_test_grid), 'x')
legend({'TX', 'RX'})
title('Test QAM-16 Constellation')


% --- Demodulate Training QAM-16 ---
rx_test_symbols = qamdemod(rx_test_grid, qam_nSymbols, UnitAveragePower=true);
if max(tx_test_symbols - rx_test_symbols) < 1e-8
    disp("Oversampled Testing receiver output matches transmitter input.");
else
    disp("Received Testing symbols do not match transmitted symbols.")
end
