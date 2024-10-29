clc; clearvars; %close all;

fcs = [19e6, 27e6, 35e6];
nFFT = 128;
M_predists = 1:7;
K_predist = 3;
upsample_amnts = [1, 2, 4, 8, 10, 15, 20, 30, 40, 60, 120];
upsample_amnt = upsample_amnts(11);

worst_train_evms = zeros(length(M_predists), 1);
worst_test_evms = zeros(length(M_predists), 1);
avg_train_evms = zeros(length(M_predists), 1);
avg_test_evms = zeros(length(M_predists), 1);

for fc=[fcs(2)]

    for ii=1:length(M_predists) 
        M_predist = M_predists(ii);
        [train_EVM_worst, test_EVM_worst, train_EVM_avg, test_EVM_avg] = get_sim( ...
            fc, nFFT, upsample_amnt, M_predist, K_predist, 'train_EVM_worst', 'test_EVM_worst', 'train_EVM_avg', 'test_EVM_avg');
        worst_train_evms(ii) = train_EVM_worst;
        worst_test_evms(ii) = test_EVM_worst;
        avg_train_evms(ii) = train_EVM_avg;
        avg_test_evms(ii) = test_EVM_avg;
        figure(1)
        plot_train_test_const(fc, nFFT, upsample_amnt, M_predist, K_predist);

        drawnow 
        pause(0.3)
    end

    m2db = @(x) 20*log10(x);
    
    figure(2)
    subplot(211)
    plot(M_predists, m2db(avg_train_evms.^2), '-o', M_predists, m2db(avg_test_evms.^2), '-*')
    title('Volt Lcl Predist. Avg. EVM')
    legend('Training', 'Test', 'Location', 'northeast')
    xlabel('M_{pred}')
    ylabel('dB')

    subplot(212)
    plot(M_predists, m2db(worst_train_evms.^2), '-o', M_predists, m2db(worst_test_evms.^2), '-*')
    title('Volt Lcl Predist. Worst EVM')
    legend('Training', 'Test', 'Location', 'northeast')
    xlabel('M_{pred}')
    ylabel('dB')
    
end


function varargout = get_sim(fc, nFFT, up_amnt, M_predist, K_predist, varargin)
    directory_name = 'D:/simulations/rls_lcl_ofdm';

     mat_filename = sprintf('rls_lcl_ofdm_fc_%g_nFFT_%g_up_%d_M_%g_K_%g.mat', fc, nFFT, up_amnt, M_predist, K_predist);
    mat_filename = strrep(mat_filename, '+', '');
    fullpath = [directory_name, '/', mat_filename];
    res = load(fullpath, varargin{:});
    varargout = struct2cell(res);
end

function plot_train_test_const(fc, nFFT, up_amnt, M_predist, K_predist)
    [tx_train_symbols, tx_test_symbols tx_train_grid, rx_train_grid, ...
        rx_test_grid, rx_train_symbols, rx_test_symbols] = get_sim( ...
            fc, nFFT, up_amnt, M_predist, K_predist, ...
            'tx_train_symbols', ...
            'tx_test_symbols', ...
            'tx_train_grid', ...
            'rx_train_grid', ...
            'rx_test_grid', ...
            'rx_train_symbols', ...
            'rx_test_symbols' ...
            );

    if max(tx_train_symbols - rx_train_symbols) < 1e-8
        disp("Oversampled Training receiver output matches transmitter input.");
        train_match = true;
    else
        disp("Received Training symbols do not match transmitted symbols.");
        train_match = false;
    end

    if max(tx_test_symbols - rx_test_symbols) < 1e-8
        disp("Oversampled Training receiver output matches transmitter input.");
        test_match = true;
    else
        disp("Received Testing symbols do not match transmitted symbols.");
        test_match = false;
    end
    
    cla
    scatter(real(tx_train_grid), imag(tx_train_grid))
    hold on
    scatter(real(rx_train_grid), imag(rx_train_grid), 'x')
    scatter(real(rx_test_grid), imag(rx_test_grid), '*')
    hold off
    
    axis(1.1.*[-1, 1, -1, 1])
    pbaspect([1, 1, 1])
    title(sprintf('Const fc=%g, nFFT=%g,  up=%g  M_{pred}=%g, K_{pred}=%g  train_mt=%d, test_mt=%d', ...
         fc, nFFT, up_amnt, M_predist, K_predist, train_match, test_match))
    box on;
    legend({'TX', 'RX_{train}', 'RX_{test}'});
end