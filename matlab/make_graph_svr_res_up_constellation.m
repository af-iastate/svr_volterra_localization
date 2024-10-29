clc; clearvars; %close all;

M_scl = 2;
fcs = [19e6, 27e6, 35e6];
nFFT = 128;
upsample_amnts = [1, 2, 4, 8, 10, 15, 20, 30, 40];
worst_train_evms = zeros(length(upsample_amnts), 1);
worst_test_evms = zeros(length(upsample_amnts), 1);
avg_train_evms = zeros(length(upsample_amnts), 1);
avg_test_evms = zeros(length(upsample_amnts), 1);

for fc=[fcs(2)]

    for ii=1:length(upsample_amnts) 
        upsample_amnt = upsample_amnts(ii);
        upfactor = 120;
        svr_factor = upfactor / upsample_amnt;
        M_svr = M_scl * (upfactor / svr_factor);
        [train_EVM_worst, test_EVM_worst, train_EVM_avg, test_EVM_avg] = get_sim( ...
            fc, nFFT, upsample_amnt, M_svr, 'train_EVM_worst', 'test_EVM_worst', 'train_EVM_avg', 'test_EVM_avg');
        worst_train_evms(ii) = train_EVM_worst;
        worst_test_evms(ii) = test_EVM_worst;
        avg_train_evms(ii) = train_EVM_avg;
        avg_test_evms(ii) = test_EVM_avg;

        figure(1)
        plot_train_test_const(fc, nFFT, upsample_amnt, M_svr)
        pause%(0.5)
    end
    
    % figure 
    % plot(upsample_amnts, avg_train_evms, '-o', upsample_amnts, avg_test_evms, '-*')
    % title(sprintf('Avg EVM: fc=%g  nFFT=%g', fc, nFFT))
    % legend('train', 'test')

end


function varargout = get_sim(fc, nFFT, up_amnt, M_svr, varargin)
    directory_name = 'simulations/svr_resolution';

    mat_filename = sprintf('svr_res_fc_%g_nFFT_%g_up_%d_Msvr_%g.mat', fc, nFFT, up_amnt, M_svr);
    mat_filename = strrep(mat_filename, '+', '');
    fullpath = [directory_name, '/', mat_filename];
    res = load(fullpath, varargin{:});
    varargout = struct2cell(res);
end

function plot_train_test_const(fc, nFFT, up_amnt, M_svr)
    [tx_train_symbols, tx_test_symbols tx_train_grid, rx_train_grid, ...
        rx_test_grid, rx_train_symbols, rx_test_symbols] = get_sim( ...
            fc, nFFT, up_amnt, M_svr, ...
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
%     title(sprintf('Const fc=%g, nFFT=%g, up=%g, M_{svr}=%g, train_mt=%d, test_mt=%d', ...
%          fc, nFFT, up_amnt, M_svr, train_match, test_match))
    title(sprintf('SVR Predist Constellation $u_d=%g$', up_amnt), 'Interpreter', 'latex')
    box on;
    set(gca, 'Position', [0.184065926831205,0.16,0.667582424817148,0.747692315795205])
    l_h = legend({'TX', 'RX_{train}', 'RX_{test}'});
    l_h.Orientation = 'horizontal';
    l_h.Position = [0.174065926831205,0.010769231640376,0.667582424817148,0.08461538374424];
    l_h.Box = false;
    l_h.EdgeColor = 'none';
end