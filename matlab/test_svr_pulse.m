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

% excitation
rng(23)
nX = 4e3;
xGauss = 1/12*randn(nX, 1);

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

% forward simulate
t = (0:nX-1)/Fs;
xPulses = [zeros(1,length(t)/2), ...
    1*sin(1e5*2*pi*t), zeros(1,length(t)/2), ...
    1*sin(2e5*2*pi*t), zeros(1,length(t)/2), ...
    1*sin(5e5*2*pi*t), zeros(1,length(t)/2), ...
    1*sin(10e5*2*pi*t), zeros(1,length(t)/2), ...
    1*sin(20e5*2*pi*t), zeros(1,length(t)/2), ...
    1*sin(50e5*2*pi*t), zeros(1,length(t)/2), ...
    ]; %+ sin(10e6*t);
xPulses = xPulses(:);

x = xPulses;
t = (0:length(x)-1)/Fs;

y = forward_f(x);

figure(1)
plot(x(5:end), y(5:end))
axis([-1, 1, -1, 1])
pbaspect([1,1,1])
title('I/O Curve (No Predistortion)')
xlabel('x(t)')
ylabel('y(t)')

figure(2)
plot(t, x)
hold on;
plot(t, y)
hold off;
title('No Predistortion')
xlabel('t')
legend({'x(t)', 'y(t)'})

% -- SVR predist --
% x_hat = volterraFilterPPM(x, h_true, M, K);
% y_hat = forward_f(x_hat);
M_svr = 17;
offset = 1;
nSamples = length(x);
Y = toeplitz(y(offset:offset+nSamples-1), [y(offset), zeros(1, M_svr-1)]);
Mdl = fitrsvm(Y, x(offset:offset+nSamples-1),  ...
...%     'Standardize', true, ...
...%       'OptimizeHyperparameters', 'auto', ...
        'KernelFunction', 'gaussian', ...
...%     'PolynomialOrder', 3, ...
    "BoxConstraint", 20, ...
    'KernelScale', 3.3132, ...
    'KKTTolerance', 5e-2, ...
    'epsilon', 1e-2);

x = xPulses;
t = (0:length(x)-1)/Fs;

X = toeplitz(x, [x(1), zeros(1, M_svr-1)]);
x_hat = Mdl.predict(X);
y_hat = forward_f(x_hat);


figure(3)
plot(x_hat(5:end), y_hat(5:end))
axis([-1, 1, -1, 1])
pbaspect([1,1,1])
title('I/O Curve (Predistortion)')
xlabel('u(t)')
ylabel('y(t)')

figure(4)
plot(t, x_hat)
hold on;
plot(t, y_hat)
hold off;
title('Predistortion')
xlabel('t')
legend({'u(t)', 'y(t)'})

for MM=0
    if MM >= 0
        xx_hat = y_hat(MM+1:end);
        xx = x(1:end-MM);
    else
        xx_hat = y_hat(1:end+MM);
        xx = x(1-MM:end);
    end

    figure(5)
    plot(t(1:end-abs(MM)), xx-xx_hat)
    title(sprintf('u(t)-x(t)'))
    xlabel('t')
    disp(norm(xx-xx_hat))
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