function [mat_filename] = simRound_multitone_thru_PA_train(M_predist, K_predist, nFFT, mat_filename)

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


% --- Init multi-tone ---
scs = 15e3;
nCycPre = 32;
upfactor = 120;
FsD = scs * nFFT/2;   % Base Sampling rate 
Fs = upfactor * FsD;
nTones = 5;
nPulses = 16;
pulseDurations = [0.01 0.01]*1e-3;

bucketSep = 1/pulseDurations(1);
numBuckets = floor(Fs/2/bucketSep);
poolSize = floor(nPulses*nTones*2.5);
freqsPool = round(bucketSep)*floor(linspace(1,numBuckets-3, poolSize));

[freqs_train, x_train, freqs_test, x_test, split_idx] = makeTonePulseDataSets(freqsPool, pulseDurations, Fs, nTones, nPulses, 10);


% --- Init Volterra RLS ---
% K_predist = 3;
% M_predist = 7; %* (upfactor / rls_factor);
lambda = 1;


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


% --- Train SVR models ---
[h, P] = volterraRlsOffline(y_train, x_train, M_predist, K_predist, 'lambda', 1);


% --- ReSimulate Training Data with Predist ---
u_train = volterraFilterPM(x_train, h, M_predist, K_predist);
y_train_hat = forward_f(u_train);


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


% --- Simulate Test Data with Predist ---
u_test = volterraFilterPM(x_test, h, M_predist, K_predist);
y_test_hat = forward_f(u_test);


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


% --- Save Results ---
save(mat_filename)


end