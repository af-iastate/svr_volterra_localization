clc; clearvars; close all;

M_predists = 1:17;
K_predist = 3;
fcs = [19e6, 27e6, 35e6];
nFFTs = [64, 128, 256];
nFFT = nFFTs(2);

worst_train_evms = zeros(length(fcs));
worst_test_evms = zeros(length(fcs));
avg_train_evms = zeros(length(fcs));
avg_test_evms = zeros(length(fcs));
train_matches = logical(zeros(length(fcs)));
test_matches = logical(zeros(length(fcs)));

for M_predist=M_predists
    for ii=1:length(fcs)
            fc_train = fcs(ii);
        
            for jj=1:length(fcs)
                fc_test = fcs(jj);
                [tx_test_grid, rx_test_grid, tx_test_symbols, rx_test_symbols] = ...
                    get_sim(fc_test, fc_train, nFFT, ...
                            M_predist, K_predist, ...
                            'tx_test_grid', ...
                            'rx_test_grid', ...
                            'tx_test_symbols', ...
                            'rx_test_symbols' ...
                    );
                test_EVMs = computeEvms(tx_test_grid, rx_test_grid);
        %         [train_EVMs, test_EVMs, train_match, test_match] = ...
        %                     plot_train_test_const(fc, nFFT, upsample_amnt, M_svr);
                worst_test_evms(ii, jj) = max(test_EVMs);
                avg_test_evms(ii, jj) = mean(test_EVMs);
                test_matches(ii, jj) = (max(tx_test_symbols - rx_test_symbols) < 1e-8);
                
            end
            
            
    end
    m2db = @(x) 20*log10(x);
    
    fprintf('M_predist = %g\n', M_predist)
%             disp(m2db(avg_test_evms))
    disp(m2db(worst_test_evms))
    disp(test_matches)
end 


function varargout = get_sim(fc_test, fc_train, nFFT, M_predist, K_predist, varargin)
    directory_name = 'simulations/rls_lcl_post_ofdm_sawp';
    mat_filename = sprintf('rls_lcl_post_ofdm_fc_test_%g_fc_train_%g_nFFT_%g_M_%g_K_%g.mat', fc_test, fc_train, nFFT, M_predist, K_predist);
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