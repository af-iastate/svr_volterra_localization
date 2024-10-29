clc; clearvars; close all;

M_scl = 2;
fcs = [19e6, 27e6, 35e6];
nFFT = 128;
upsample_amnts = [1, 2, 4, 8, 10, 15, 20, 30, 40];
worst_train_evms = zeros(length(upsample_amnts), 1);
worst_test_evms = zeros(length(upsample_amnts), 1);
avg_train_evms = zeros(length(upsample_amnts), 1);
avg_test_evms = zeros(length(upsample_amnts), 1);
train_matches = logical(zeros(length(upsample_amnts), 1));
test_matches = logical(zeros(length(upsample_amnts), 1));

for fc=[fcs(2)]
    for ii=1:length(upsample_amnts)

        upsample_amnt = upsample_amnts(ii);
        upfactor = 120;
        svr_factor = upfactor / upsample_amnt;
        M_svr = M_scl * (upfactor / svr_factor);

        figure(1)
        [train_EVMs, test_EVMs, train_match, test_match] = plot_train_test_const(fc, nFFT, upsample_amnt, M_svr);
        worst_train_evms(ii) = max(train_EVMs);
        worst_test_evms(ii) = max(test_EVMs);
        avg_train_evms(ii) = mean(train_EVMs);
        avg_test_evms(ii) = mean(test_EVMs);
        train_matches(ii) = logical(train_match);
        test_matches(ii) = logical(test_match);
        
        drawnow
        pause(0.3)
    end
    
    m2db = @(x) 20*log10(x);

    figure(2)
    C = colororder;

    subplot(211)
    cla
    hold on

    % --- Train ---
    h_train = plot(upsample_amnts, m2db(avg_train_evms), 'color', C(1,:));
    scatter(upsample_amnts(~train_matches), m2db(avg_train_evms(~train_matches)), 'CData', h_train.Color, ...
        'SizeData', 25, 'MarkerFaceColor', 'White');
    scatter(upsample_amnts(train_matches), m2db(avg_train_evms(train_matches)), 'Cdata', h_train.Color, ...
        'marker', 'o', 'MarkerFaceColor', h_train.Color, 'MarkerEdgeColor', 'none', 'SizeData', 25)

    % --- Test ---
    h_test = plot(upsample_amnts, m2db(avg_test_evms), 'color', C(2,:));
    scatter(upsample_amnts(~test_matches), m2db(avg_test_evms(~test_matches)), 'CData', h_test.Color, ...
        'SizeData', 25, 'MarkerFaceColor', 'White');
    scatter(upsample_amnts(test_matches), m2db(avg_test_evms(test_matches)), 'Cdata', h_test.Color, ...
        'marker', 'o', 'MarkerFaceColor', h_test.Color, 'MarkerEdgeColor', 'none', 'SizeData', 25)
    
    title('SVR Predist. Avg. EVM', 'Interpreter', 'latex')
    legend({'Training', '', '', 'Test'}, 'Location', 'north', 'Orientation', 'horizontal', 'Interpreter', 'latex')
    xlabel('$u_d$', 'Interpreter', 'latex')
    ylabel('dB', 'Interpreter', 'latex')

    hold off


    subplot(212)
    cla
    hold on

    % --- Train ---
    h_train = plot(upsample_amnts, m2db(worst_train_evms), 'color', C(1,:));
    scatter(upsample_amnts(~train_matches), m2db(worst_train_evms(~train_matches)), 'CData', h_train.Color, ...
        'SizeData', 25, 'MarkerFaceColor', 'White');
    scatter(upsample_amnts(train_matches), m2db(worst_train_evms(train_matches)), 'Cdata', h_train.Color, ...
        'marker', 'o', 'MarkerFaceColor', h_train.Color, 'MarkerEdgeColor', 'none', 'SizeData', 25)

    % --- Test ---
    h_test = plot(upsample_amnts, m2db(worst_test_evms), 'color', C(2,:));
    scatter(upsample_amnts(~test_matches), m2db(worst_test_evms(~test_matches)), 'CData', h_test.Color, ...
        'SizeData', 25, 'MarkerFaceColor', 'White');
    scatter(upsample_amnts(test_matches), m2db(worst_test_evms(test_matches)), 'Cdata', h_test.Color, ...
        'marker', 'o', 'MarkerFaceColor', h_test.Color, 'MarkerEdgeColor', 'none', 'SizeData', 25)
    
    
    title('SVR Predist. Worst EVM', 'Interpreter', 'latex')
    legend({'Training', '', '', 'Test'}, 'Location', 'north', 'Orientation', 'horizontal', 'Interpreter', 'latex')
    xlabel('$u_d$', 'Interpreter', 'latex')
    ylabel('dB', 'Interpreter', 'latex')
end


function varargout = get_sim(fc, nFFT, up_amnt, M_svr, varargin)
    directory_name = 'D:/simulations/svr_resolution';

    mat_filename = sprintf('svr_res_fc_%g_nFFT_%g_up_%d_Msvr_%g.mat', fc, nFFT, up_amnt, M_svr);
    mat_filename = strrep(mat_filename, '+', '');
    fullpath = [directory_name, '/', mat_filename];
    res = load(fullpath, varargin{:});
    varargout = struct2cell(res);
end

function [train_EVMs, test_EVMs, train_match, test_match] = plot_train_test_const(fc, nFFT, up_amnt, M_svr)
    [tx_train_symbols, tx_test_symbols, tx_train_grid, tx_test_grid, rx_train_grid, ...
        rx_test_grid, rx_train_symbols, rx_test_symbols] = get_sim( ...
            fc, nFFT, up_amnt, M_svr, ...
            'tx_train_symbols', ...
            'tx_test_symbols', ...
            'tx_train_grid', ...
            'tx_test_grid', ...
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
    
    train_EVMs = computeEvms(tx_train_grid, rx_train_grid);
    test_EVMs = computeEvms(tx_test_grid, rx_test_grid);

    cla
    scatter(real(tx_train_grid), imag(tx_train_grid))
    hold on
    scatter(real(rx_train_grid), imag(rx_train_grid), 'x')
    scatter(real(rx_test_grid), imag(rx_test_grid), '*')
    hold off
    
    axis(1.1.*[-1, 1, -1, 1])
    pbaspect([1, 1, 1])
    title(sprintf('Const fc=%g, nFFT=%g, up=%g, M_{svr}=%g, train_mt=%d, test_mt=%d', ...
         fc, nFFT, up_amnt, M_svr, train_match, test_match))
    box on;
    legend({'TX', 'RX_{train}', 'RX_{test}'});
end


function [evms] = computeEvms(tx_grid, rx_grid)
    evms = abs(tx_grid - rx_grid) ./ abs(tx_grid);
end