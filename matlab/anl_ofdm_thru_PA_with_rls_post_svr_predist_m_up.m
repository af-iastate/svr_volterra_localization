clc; clearvars; close all;

fcs = [19e6, 27e6, 35e6];
nFFTs = [64, 128, 256];
nFFT = nFFTs(2);
M_svrs = 1:17;

worst_train_evms = zeros(length(M_svrs), 1);
worst_test_evms = zeros(length(M_svrs), 1);
avg_train_evms = zeros(length(M_svrs), 1);
avg_test_evms = zeros(length(M_svrs), 1);
train_matches = logical(zeros(length(M_svrs), 1));
test_matches = logical(zeros(length(M_svrs), 1));

figure(1)
clf

for jj=1:length(fcs)
    fc = fcs(jj);

    for ii=1:length(M_svrs) 
        M_svr = M_svrs(ii);
    
        figure(1)
        [train_EVMs, test_EVMs, train_match, test_match] = plot_train_test_const(fc, nFFT, M_svr);

        worst_train_evms(ii) = max(train_EVMs);
        worst_test_evms(ii) = max(test_EVMs);
        avg_train_evms(ii) = mean(train_EVMs);
        avg_test_evms(ii) = mean(test_EVMs);
        train_matches(ii) = logical(train_match);
        test_matches(ii) = logical(test_match);

        drawnow 
    end
    
    m2db = @(x) 20*log10(x);

    C = colororder;
    
    figure(1)
    subplot(3, 1, jj)
    cla;
    hold on
    
    % --- Train ---
    h_train = plot(M_svrs, m2db(worst_train_evms), "Color", C(1,:));
    scatter(M_svrs(~train_matches), m2db(worst_train_evms(~train_matches)), 'CData', h_train.Color, ...
        'SizeData', 25, 'MarkerFaceColor', 'White');
    scatter(M_svrs(train_matches), m2db(worst_train_evms(train_matches)), 'Cdata', h_train.Color, ...
        'marker', 'o', 'MarkerFaceColor', h_train.Color, 'MarkerEdgeColor', 'none', 'SizeData', 25)
    % --- Test ---
    h_test = plot(M_svrs, m2db(worst_test_evms), "Color", C(2,:));
    scatter(M_svrs(~test_matches), m2db(worst_test_evms(~test_matches)), 'CData', h_test.Color, ...
        'SizeData', 25, 'MarkerFaceColor', 'White');
    scatter(M_svrs(test_matches), m2db(worst_test_evms(test_matches)), 'Cdata', h_test.Color, ...
        'marker', 'o', 'MarkerFaceColor', h_test.Color, 'MarkerEdgeColor', 'none', 'SizeData', 25)
    hold off
    my_loc = 'north';
    if any([m2db(worst_train_evms(6:12)); m2db(worst_test_evms(6:12))] >= -6)
        my_loc = 'south';
    end
    title(sprintf('Worst EVM $f_c=%d$\\,MHz', fc/1e6), 'Interpreter','latex')
    xlabel('$M_\textrm{pred}$', 'Interpreter', 'latex')
    ylabel('dB', 'Interpreter','latex')
    xticks(M_svrs)
    axis([xlim, -44, 0])
    yticks(-40:10:0)
    legend({'Train', '', '', 'Test'}, 'Orientation', 'horizontal', 'Location',my_loc)

end
drawnow


function varargout = get_sim(fc, nFFT, M_svr, varargin)
    directory_name = 'simulations/svr_lcl_post_ofdm';
    mat_filename = sprintf('svr_lcl_post_ofdm_fc_%g_nFFT_%g_Msvr_%g.mat', fc, nFFT, M_svr);
    mat_filename = strrep(mat_filename, '+', '');
    mat_filename = strrep(mat_filename, '+', '');
    fullpath = [directory_name, '/', mat_filename];
    res = load(fullpath, varargin{:});
    varargout = struct2cell(res);
end


function [train_EVMs, test_EVMs, train_match, test_match] = plot_train_test_const(fc, nFFT, M_svr)
    [tx_train_symbols, tx_test_symbols, tx_train_grid, tx_test_grid, rx_train_grid, ...
        rx_test_grid, rx_train_symbols, rx_test_symbols] = get_sim( ...
            fc, nFFT, M_svr, ...
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
end



function [evms] = computeEvms(tx_grid, rx_grid)
    evms = abs(tx_grid - rx_grid) ./ abs(tx_grid);
end