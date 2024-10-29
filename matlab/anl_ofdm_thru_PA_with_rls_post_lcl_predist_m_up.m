clc; clearvars; %close all;

fcs = [19e6, 27e6, 35e6];
nFFTs = [64, 128, 256];
nFFT = nFFTs(2);
M_predists = 1:17;
K_predist = 3;

worst_train_evms = zeros(length(M_predists), 1);
worst_test_evms = zeros(length(M_predists), 1);
avg_train_evms = zeros(length(M_predists), 1);
avg_test_evms = zeros(length(M_predists), 1);
train_matches = logical(zeros(length(M_predists), 1));
test_matches = logical(zeros(length(M_predists), 1));

figure(1)
clf
for jj=1:length(fcs)
    fc = fcs(jj);

    for ii=1:length(M_predists) 
        M_predist = M_predists(ii);
    
        figure(1)
        [train_EVMs, test_EVMs, train_match, test_match] = plot_train_test_const(fc, nFFT, M_predist, K_predist);

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

    C = colororder;
    
    figure(1)
    subplot(3, 1, jj)
    cla;
    hold on
    
    % --- Train ---
    h_train = plot(M_predists, m2db(worst_train_evms), "Color", C(1,:));
    scatter(M_predists(~train_matches), m2db(worst_train_evms(~train_matches)), 'CData', h_train.Color, ...
        'SizeData', 25, 'MarkerFaceColor', 'White');
    scatter(M_predists(train_matches), m2db(worst_train_evms(train_matches)), 'Cdata', h_train.Color, ...
        'marker', 'o', 'MarkerFaceColor', h_train.Color, 'MarkerEdgeColor', 'none', 'SizeData', 25)
    % --- Test ---
    h_test = plot(M_predists, m2db(worst_test_evms), "Color", C(2,:));
    scatter(M_predists(~test_matches), m2db(worst_test_evms(~test_matches)), 'CData', h_test.Color, ...
        'SizeData', 25, 'MarkerFaceColor', 'White');
    scatter(M_predists(test_matches), m2db(worst_test_evms(test_matches)), 'Cdata', h_test.Color, ...
        'marker', 'o', 'MarkerFaceColor', h_test.Color, 'MarkerEdgeColor', 'none', 'SizeData', 25)
    hold off
    my_loc = 'north';
    if any([m2db(worst_train_evms(6:12)); m2db(worst_test_evms(6:12))] >= -6)
        my_loc = 'south';
    end
    title(sprintf('Worst EVM $f_c=%d$\\,MHz', fc/1e6), 'Interpreter','latex')
    xlabel('$M_\textrm{pred}$', 'Interpreter', 'latex')
    ylabel('dB', 'Interpreter','latex')
    xticks(M_predists)
    axis([xlim, -44, 0])
    yticks(-40:10:0)
    legend({'Train', '', '', 'Test'}, 'Orientation', 'horizontal', 'Location',my_loc)
    
end

% gif_path = sprintf('gif/rls_post_lcl_predist_m_up/evms.gif');
% exportgraphics(gcf, gif_path, "Append", true)


function varargout = get_sim(fc, nFFT, M_predist, K_predist, varargin)
    directory_name = 'simulations/rls_lcl_post_ofdm';
    mat_filename = sprintf('rls_lcl_post_ofdm_fc_%g_nFFT_%g_M_%g_K_%g.mat', fc, nFFT, M_predist, K_predist);
    mat_filename = strrep(mat_filename, '+', '');
    fullpath = [directory_name, '/', mat_filename];
    res = load(fullpath, varargin{:});
    varargout = struct2cell(res);
end


function [train_EVMs, test_EVMs, train_match, test_match] = plot_train_test_const(fc, nFFT, M_predist, K_predist)
    [tx_train_symbols, tx_test_symbols, tx_train_grid, tx_test_grid, rx_train_grid, ...
        rx_test_grid, rx_train_symbols, rx_test_symbols] = get_sim( ...
            fc, nFFT, M_predist, K_predist, ...
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
    
%     figure(2)
%     cla
%     scatter(real(tx_train_grid), imag(tx_train_grid))
%     hold on
%     scatter(real(rx_train_grid), imag(rx_train_grid), 'x')
%     scatter(real(rx_test_grid), imag(rx_test_grid), 'x')
%     hold off
%     
%     axis(1.1.*[-1, 1, -1, 1])
%     pbaspect([1, 1, 1])
%     title(sprintf('Const fc=%d MHz, nFFT=%g, M_{pred}=%02g, K_{pred}=%g  train_{mt}=%d, test_{mt}=%d', ...
%          round(fc/1e6), nFFT, M_predist, K_predist, train_match, test_match))
%     box on;
%     legend({'TX', 'RX_{train}', 'RX_{test}'});
    
%     gif_path = sprintf('gif/rls_post_lcl_predist_m_up/fc_%d.gif', round(fc/1e6));
%     exportgraphics(gca, gif_path, "Append", true)
end



function [evms] = computeEvms(tx_grid, rx_grid)
    evms = abs(tx_grid - rx_grid) ./ abs(tx_grid);
end