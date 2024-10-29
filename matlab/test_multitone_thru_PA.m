clearvars; close all; clc;
% load Volterra filtering method scripts
currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
addpath(fullfile(pathstr, 'ofdm'));
addpath(fullfile(pathstr, 'known-inverse'));
addpath(fullfile(pathstr, 'methods/morhac'));
addpath(fullfile(pathstr, 'methods/pm'));
addpath(fullfile(pathstr, 'methods/ppm'));


% --- Init multi-tone ---
Fs = 115.2e6;
nTones = 5;
nPulses = 16;
pulseDurations = [0.01 0.01]*1e-3;

bucketSep = 1/pulseDurations(1);
numBuckets = floor(Fs/2/bucketSep);
poolSize = floor(nPulses*nTones*2.5);
freqsPool = round(bucketSep)*floor(linspace(1,numBuckets-3, poolSize));

[freqs_train, x_train, freqs_test, x_test, split_idx] = makeTonePulseDataSets(freqsPool, pulseDurations, Fs, nTones, nPulses, 10);


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


% --- Plot OFDM RS/TX signal ---
figure(1)
t_x_train = (0:length(x_train)-1)./Fs;
t_y_train = (0:length(y_train)-1)./Fs;

subplot(211)
plot(t_x_train, x_train)
hold on
plot(t_y_train, y_train)
legend({'TX', 'RX'})
title('RX/TX Time domain')

subplot(212)
pspectrum(x_train, Fs)
hold on
pspectrum(y_train, Fs)
legend({'TX', 'RX'})