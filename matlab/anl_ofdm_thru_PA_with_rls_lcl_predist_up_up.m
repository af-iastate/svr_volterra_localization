clc; clearvars; %close all;

fcs = [19e6, 27e6, 35e6];
nFFTs = [64, 128, 256];
nFFT = nFFTs(2);
M_predists = 1:8;
M_predist = 1;
K_predist = 3;
upsample_amnts = [1, 2, 4, 8, 10, 15, 20, 30, 40, 60, 120];
upsample_amnt = upsample_amnts(9);

worst_train_evms = zeros(length(M_predists), 1);
worst_test_evms = zeros(length(M_predists), 1);
avg_train_evms = zeros(length(M_predists), 1);
avg_test_evms = zeros(length(M_predists), 1);
train_matches = logical(zeros(length(M_predists), 1));
test_matches = logical(zeros(length(M_predists), 1));



figure(1)
clf
for M_predist=M_predists(2)
    for jj=1:length(fcs)
        fc = fcs(jj);
    
        for ii=1:length(upsample_amnts) 
            upsample_amnt = upsample_amnts(ii);
        
            figure(1)
            [train_EVMs, test_EVMs, train_match, test_match] = plot_train_test_const(fc, nFFT, upsample_amnt, M_predist, K_predist);
    
            worst_train_evms(ii) = max(train_EVMs);
            worst_test_evms(ii) = max(test_EVMs);
            avg_train_evms(ii) = mean(train_EVMs);
            avg_test_evms(ii) = mean(test_EVMs);
            train_matches(ii) = logical(train_match);
            test_matches(ii) = logical(test_match);
    
            drawnow 
    %         pause(0.3)
        end
        
        m2db = @(x) 20*log10(x);
        tic_locs = 1:length(upsample_amnts);
        C = colororder;
    
        figure(1)
        subplot(3, 1, jj)
        cla;
        hold on
        % --- Train ---
        h_train = plot(tic_locs, m2db(worst_train_evms), 'Color', C(1,:));
        scatter(tic_locs(~train_matches), m2db(worst_train_evms(~train_matches)), 'CData', h_train.Color, ...
            'SizeData', 25, 'MarkerFaceColor', 'White');
        scatter(tic_locs(train_matches), m2db(worst_train_evms(train_matches)), 'Cdata', h_train.Color, ...
            'marker', 'o', 'MarkerFaceColor', h_train.Color, 'MarkerEdgeColor', 'none', 'SizeData', 25)
        % --- Test ---
        h_test = plot(tic_locs, m2db(worst_test_evms), 'Color', C(2,:));
        scatter(tic_locs(~test_matches), m2db(worst_test_evms(~test_matches)), 'CData', h_test.Color, ...
            'SizeData', 25, 'MarkerFaceColor', 'White');
        scatter(tic_locs(test_matches), m2db(worst_test_evms(test_matches)), 'Cdata', h_test.Color, ...
            'marker', 'o', 'MarkerFaceColor', h_test.Color, 'MarkerEdgeColor', 'none', 'SizeData', 25)
        hold off
        title(sprintf('Worst EVM $M_\\mathrm{pred}=%d$, $f_c=%d$\\,MHz', M_predist, fc/1e6), 'Interpreter','latex')
        xlabel('$u_d$', 'Interpreter', 'latex')
        ylabel('dB', 'Interpreter','latex')
        xticks(tic_locs)
        xticklabels(upsample_amnts)
        axis([xlim, -35, 10])
        legend({'Train', '', '', 'Test'}, 'Orientation', 'horizontal', 'Location','north')
    end
    drawnow
    pause
    
%     gif_path = sprintf('gif/rls_lcl_predist_up_up/evms.gif');
%     exportgraphics(gcf, gif_path, "Append", true)
end

function varargout = get_sim(fc, nFFT, upsample_amnt, M_predist, K_predist, varargin)
    directory_name = 'simulations/rls_lcl_ofdm';
    mat_filename = sprintf('rls_lcl_ofdm_fc_%g_nFFT_%g_up_%d_M_%g_K_%g.mat', fc, nFFT, upsample_amnt, M_predist, K_predist);
    mat_filename = strrep(mat_filename, '+', '');
    fullpath = [directory_name, '/', mat_filename];
    res = load(fullpath, varargin{:});
    varargout = struct2cell(res);
end


function [train_EVMs, test_EVMs, train_match, test_match] = plot_train_test_const(fc, nFFT, upsample_amnt, M_predist, K_predist)
    [tx_train_symbols, tx_test_symbols, tx_train_grid, tx_test_grid, rx_train_grid, ...
        rx_test_grid, rx_train_symbols, rx_test_symbols] = get_sim( ...
            fc, nFFT, upsample_amnt, M_predist, K_predist, ...
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
%     title(sprintf('Const fc=%d MHz, nFFT=%g, up=%03d, M_{pred}=%g, K_{pred}=%g  train_{mt}=%d, test_{mt}=%d', ...
%          round(fc/1e6), nFFT, upsample_amnt, M_predist, K_predist, train_match, test_match))
%     box on;
%     legend({'TX', 'RX_{train}', 'RX_{test}'});
    
%     gif_path = sprintf('gif/rls_lcl_predist_up_up/fc_%d_Mpred_%d.gif', round(fc/1e6), M_predist);
%     exportgraphics(gca, gif_path, "Append", true)
end



function [evms] = computeEvms(tx_grid, rx_grid)
    evms = abs(tx_grid - rx_grid) ./ abs(tx_grid);
end