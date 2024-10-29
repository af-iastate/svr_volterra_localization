function [mat_filename] = simRound_swap_ofdm_thru_PA_with_rls_post_lcl_predist(fc_test, fc_train, nFFT, M_predist, K_predist, mat_filename)

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
% fc = 20e6;
scs = 15e3;
% nFFT = 2^7;
nCycPre = 32;
upfactor = 120;
qam_bitsPerSymbol = 4;               % Bits per symbol
qam_nSymbols = 2^qam_bitsPerSymbol;  % QAM-16
headroom_coeff = 1.5;


% --- Init Volterra RLS ---
% K_predist = 3;
% M_predist = 7; %* (upfactor / rls_factor);
lambda = 1;


% --- Init Training ---
directory_name = 'D:/simulations/rls_lcl_post_ofdm';
sim_mat_filename = sprintf('rls_lcl_post_ofdm_fc_%g_nFFT_%g_M_%g_K_%g.mat', fc_train, nFFT, M_predist, K_predist);
sim_mat_filename = strrep(sim_mat_filename, '+', '');
fullpath = [directory_name, '/', sim_mat_filename];
res = load(fullpath, 'h');
h = res.h;
clearvars res;


% --- Init Testing ---
rng(29)
tx_test_symbols = [0, 0, randi([0, qam_nSymbols-1], nFFT - 4, 1).', 0, 0].';
tx_test_grid = qammod(tx_test_symbols, qam_nSymbols, UnitAveragePower=true);

[x_test_iq, Fs] = ofdm_modulate(tx_test_grid, scs, 'nCycPre', nCycPre, 'upfactor', upfactor);
x_test = iq_modulate(x_test_iq, Fs, fc_test);


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


% --- Simulate Test Data with Predist ---
u_test = volterraFilterPM(x_test, h, M_predist, K_predist);
y_test_hat = forward_f(u_test);
y_test_hat_iq = iq_demodulate(y_test_hat, Fs, fc_test, (headroom_coeff*scs*nFFT/2));


% --- Plot Training OFDM RS/TX with Predistortion ---
figure
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
figure
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
figure
scatter(real(tx_test_grid), imag(tx_test_grid))
hold on
scatter(real(rx_test_grid), imag(rx_test_grid), 'x')
legend({'TX', 'RX'})
title('Test QAM-16 Constellation')

test_EVMs = sqrt(abs(tx_test_grid - rx_test_grid));
test_EVM_avg = sum(test_EVMs) / length(tx_test_grid);
test_EVM_worst = max(test_EVMs);

% --- Demodulate Training QAM-16 ---
rx_test_symbols = qamdemod(rx_test_grid, qam_nSymbols, UnitAveragePower=true);
if max(tx_test_symbols - rx_test_symbols) < 1e-8
    disp("Oversampled Testing receiver output matches transmitter input.");
else
    disp("Received Testing symbols do not match transmitted symbols.")
end

fprintf('Test EVM avg = %g\n', test_EVM_avg)
fprintf('Test EVM wrst = %g\n\n', test_EVM_worst)

% --- Save Results ---
save(mat_filename)


end