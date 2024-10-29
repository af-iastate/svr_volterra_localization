clearvars; close all; clc;

% load Volterra filtering method scripts
currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
addpath(fullfile(pathstr, 'known-inverse'));
addpath(fullfile(pathstr, 'methods/morhac'));
addpath(fullfile(pathstr, 'methods/pm'));
addpath(fullfile(pathstr, 'methods/ppm'));

len = @length;
Fs = 125e6;

%% --- Excitation Corner ---
rng(23)

% Gaussian
nX = 4e3;
xGauss = 1/12*randn(nX, 1);

% Single tone
t_training = (0:nX-1)/Fs;
xPulses = [zeros(1,length(t_training)/2), ...
    1*sin(1e5*2*pi*t_training), zeros(1,length(t_training)/2), ...
    1*sin(2e5*2*pi*t_training), zeros(1,length(t_training)/2), ...
    1*sin(5e5*2*pi*t_training), zeros(1,length(t_training)/2), ...
    1*sin(10e5*2*pi*t_training), zeros(1,length(t_training)/2), ...
    1*sin(20e5*2*pi*t_training), zeros(1,length(t_training)/2), ...
    1*sin(50e5*2*pi*t_training), zeros(1,length(t_training)/2), ...
    ]; %+ sin(10e6*t);
xPulses = xPulses(:);

% 5-tone pulse
durationTwoToneAdditional = [0.01 0.01]*1e-3;
bucketSep = 1/durationTwoToneAdditional(1);
numBuckets = floor(Fs/2/bucketSep);
poolSize = floor(16*5*2.5);
freqsPool = round(bucketSep)*floor(linspace(1,numBuckets-3, poolSize));
[freqsTwoTone, xMultiTone, freqsTwoToneAlt, xMultiToneAlt, splitIdx] = ...
    makeTonePulseDataSets(freqsPool, durationTwoToneAdditional, Fs, 5, 16);


%% --- System Definitions ---
memblock_a = 0.999*GenBandstopIr(10, Fs, 30e6, 40e6);
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


%% --- Perform Simulations ---
x_training = xMultiTone;
t_training = (0:length(x_training)-1)/Fs;

y_training = forward_f(x_training);

figure(1)
plot(x_training(5:end), y_training(5:end))
axis([-1, 1, -1, 1])
pbaspect([1,1,1])
title('I/O Curve (No Predistortion)')
xlabel('x(t)')
ylabel('y(t)')

figure(2)
plot(t_training, x_training)
hold on;
plot(t_training, y_training)
hold off;
title('No Predistortion')
xlabel('t')
legend({'x(t)', 'y(t)'})

%% --- SVR Train ---

% x_hat = volterraFilterPPM(x, h_true, M, K);
% y_hat = forward_f(x_hat);

M_svr = 30;
offset = 1;
nSamples = length(x_training);
Y_training = toeplitz(y_training(offset:offset+nSamples-1), [y_training(offset), zeros(1, M_svr-1)]);

Mdl = fitrsvm(Y_training, x_training(offset:offset+nSamples-1),  ...
...%     'Standardize', true, ...
...%      'OptimizeHyperparameters', 'auto', ...
     'KernelFunction', 'gaussian', ...
...%      'PolynomialOrder', 3, ...
     "BoxConstraint", 7, ...
     'KernelScale', 0.5*M_svr, ...
    'KKTTolerance', 7e-3, ...
    'epsilon', 1e-3);

X_training = toeplitz(x_training, [x_training(1), zeros(1, M_svr-1)]);
u_training = Mdl.predict(X_training);
y_hat = forward_f(u_training);

figure(3)
plot(u_training(5:end), y_hat(5:end))
axis([-1, 1, -1, 1])
pbaspect([1,1,1])
title('I/O Curve (Predistortion Training)')
xlabel('u(t)')
ylabel('y(t)')

figure(4)
plot(t_training, x_training)
hold on;
plot(t_training, y_hat)
hold off;
title('Predistortion Training')
xlabel('t')
legend({'x(t)', 'y(t)'})

for MM=0
    if MM >= 0
        yy_hat = y_hat(MM+1:end);
        xx_training = x_training(1:end-MM);
    else
        yy_hat = y_hat(1:end+MM);
        xx_training = x_training(1-MM:end);
    end

    figure(5)
    plot(t_training(1:end-abs(MM)), xx_training-yy_hat)
    title(sprintf('x(t)-y(t) Training'))
    xlabel('t')
    disp(norm(xx_training-yy_hat))
    pause(0.5)
end


%% --- SVR Test ---
x_test = xMultiToneAlt;
t_test = (0:length(x_test)-1)/Fs;

X_test = toeplitz(x_test, [x_test(1), zeros(1, M_svr-1)]);
u_test = Mdl.predict(X_test);
y_hat = forward_f(u_test);

figure(6)
plot(u_test(5:end), y_hat(5:end))
axis([-1, 1, -1, 1])
pbaspect([1,1,1])
title('I/O Curve (Predistortion Test)')
xlabel('u(t)')
ylabel('y(t)')

figure(7)
plot(t_test, x_test)
hold on;
plot(t_test, y_hat)
hold off;
title('Predistortion Test')
xlabel('t')
legend({'x(t)', 'y(t)'})

for MM=0
    if MM >= 0
        yy_hat = y_hat(MM+1:end);
        xx_test = x_test(1:end-MM);
    else
        yy_hat = y_hat(1:end+MM);
        xx_test = x_test(1-MM:end);
    end

    figure(8)
    plot(t_test(1:end-abs(MM)), xx_test-yy_hat)
    title(sprintf('x(t)-y(t) Test'))
    xlabel('t')
    disp(norm(xx_test-yy_hat))
    pause(0.5)
end



%X = X(M:end, :);

% Y = getXMatrix(y, M, K);
% h_hat = pinv(Y) * x;
% 
% figure(6)
% stem(h_true)
% hold on
% stem(h_hat)
% 
% figure(7)
% stem(h_true-h_hat)
% title(sprintf('Error, max err=%g', norm(h_true-h_hat)))