function [mat_filename] = simRound_ofdm_thru_PA_with_svr_predist(fc, upsample_amnt, nFFT, M_svr, mat_filename)
close all;

% load Volterra filtering method scripts
currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
addpath(fullfile(pathstr, 'ofdm'));
addpath(fullfile(pathstr, 'known-inverse'));
addpath(fullfile(pathstr, 'methods/morhac'));
addpath(fullfile(pathstr, 'methods/pm'));
addpath(fullfile(pathstr, 'methods/ppm'));


% --- Init OFDM ---
% fc = 20e6;
scs = 15e3;
% nFFT = 128;
nCycPre = 32;
upfactor = 120;
qam_bitsPerSymbol = 4;               % Bits per symbol
qam_nSymbols = 2^qam_bitsPerSymbol;  % QAM-16
headroom_coeff = 1.5;


% --- Init SVR ---
svr_factor = upfactor / upsample_amnt;
% M_svr = M_scl * (upfactor / svr_factor);
svr_options = { ...
%     'Standardize', true
    'OptimizeHyperparameters', 'auto'
    'KernelFunction', 'gaussian' 
%     'PolynomialOrder', 3
%     'BoxConstraint', 1 
%     'KernelScale', 1*M_svr
    'KKTTolerance', 7e-4
    'epsilon', 5e-5 
}.';

f_d = @(x) downsample(x, svr_factor);
f_u = @(x) sinc_upsample(x, svr_factor);


% --- Init Training ---
rng(29)% + fc + upsample_amnt + 10* round(sqrt(nFFT)) + M_scl*1024)
tx_train_symbols = [0, 0, randi([0, qam_nSymbols-1], nFFT - 4, 1).', 0, 0].';
tx_train_grid = qammod(tx_train_symbols, qam_nSymbols, UnitAveragePower=true);

[x_train_iq, Fs] = ofdm_modulate(tx_train_grid, scs, 'nCycPre', nCycPre, 'upfactor', upfactor);
x_train_iq_ds = f_d(x_train_iq);
x_train = iq_modulate(x_train_iq, Fs, fc);


% --- Init Testing ---
rng(188)% + fc + upsample_amnt + 10* round(sqrt(nFFT)) + M_scl*1024)
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
figure
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
figure
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
offset = 1 + nCycPre*(upfactor / svr_factor);
nSamples = nFFT * (upfactor / svr_factor);
Y_train_iq_ds = toeplitz(y_train_iq_ds(offset:offset+nSamples-1), [y_train_iq_ds(offset), zeros(1, M_svr-1)]);

I_Mdl = fitrsvm([real(Y_train_iq_ds), imag(Y_train_iq_ds)], real(x_train_iq_ds(offset:offset+nSamples-1)), svr_options{:});
h_I_Mdl = gcf();
h_I_Mdl_struct = handle2struct(h_I_Mdl);
clearvars h_I_Mdl;

Q_Mdl = fitrsvm([real(Y_train_iq_ds), imag(Y_train_iq_ds)], imag(x_train_iq_ds(offset:offset+nSamples-1)), svr_options{:});
h_Q_Mdl = gcf();
h_Q_Mdl_struct = handle2struct(h_Q_Mdl);
clearvars h_Q_Mdl;

% --- ReSimulate Training Data with Predist ---
X_train_iq_ds = toeplitz(x_train_iq_ds, [x_train_iq_ds(1), zeros(1, M_svr-1)]);
u_train_iq = f_u(I_Mdl.predict([real(X_train_iq_ds), imag(X_train_iq_ds)])) ...
    + 1j*f_u(Q_Mdl.predict([real(X_train_iq_ds), imag(X_train_iq_ds)]));
u_train = iq_modulate(u_train_iq, Fs, fc);
y_train_hat = forward_f(u_train);
y_train_hat_iq = iq_demodulate(y_train_hat, Fs, fc, (headroom_coeff*scs*nFFT/2));


% --- Plot Training OFDM RS/TX with Predistortion ---
figure
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
figure
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
figure
scatter(real(tx_train_grid), imag(tx_train_grid))
hold on
scatter(real(rx_train_grid), imag(rx_train_grid), 'x')
legend({'TX', 'RX'})
title('Training QAM-16 Constellation')

train_EVMs = sqrt(abs(tx_train_grid - rx_train_grid));
train_EVM_avg = sum(train_EVMs) / length(tx_train_grid);
train_EVM_worst = max(train_EVMs);

% --- Demodulate Training QAM-16 ---
rx_train_symbols = qamdemod(rx_train_grid, qam_nSymbols, UnitAveragePower=true);
if max(tx_train_symbols - rx_train_symbols) < 1e-8
    disp("Oversampled Training receiver output matches transmitter input.");
else
    disp("Received Training symbols do not match transmitted symbols.")
end

fprintf('Train EVM avg = %g\n', train_EVM_avg)
fprintf('Train EVM wrst = %g\n\n', train_EVM_worst)


% --- Simulate Test Data with Predist ---
X_test_iq_ds = toeplitz(x_test_iq_ds, [x_test_iq_ds(1), zeros(1, M_svr-1)]);
u_test_iq = f_u(I_Mdl.predict([real(X_test_iq_ds), imag(X_test_iq_ds)])) ...
    + 1j*f_u(Q_Mdl.predict([real(X_test_iq_ds), imag(X_test_iq_ds)]));
u_test = iq_modulate(u_test_iq, Fs, fc);
y_test_hat = forward_f(u_test);
y_test_hat_iq = iq_demodulate(y_test_hat, Fs, fc, (headroom_coeff*scs*nFFT/2));


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