clc; clearvars; close all;

M_svrs = 1:17;
fcs = [19e6, 27e6, 35e6];
nFFTs = [64, 128, 256];
nFFT = nFFTs(2);

worst_train_evms = zeros(length(fcs));
worst_test_evms = zeros(length(fcs));
avg_train_evms = zeros(length(fcs));
avg_test_evms = zeros(length(fcs));
train_matches = logical(zeros(length(fcs)));
test_matches = logical(zeros(length(fcs)));

for M_svr=M_svrs
    for ii=1:length(fcs)
            fc_train = fcs(ii);
        
            for jj=1:length(fcs)
                fc_test = fcs(jj);
                [tx_test_grid, rx_test_grid, tx_test_symbols, rx_test_symbols] = ...
                    get_sim(fc_test, fc_train, nFFT, ...
                            M_svr, ...
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
    
    fprintf('M_pred = %g\n', M_svr)
%             disp(m2db(avg_test_evms))
    disp(m2db(worst_test_evms))
    disp(test_matches)

    figure(1)
    clf
    [X,Y] = meshgrid(1:3, 1:3);
    Z = m2db(worst_test_evms);
    b_h = bar3(Z);
    axis([0, 4, 0, 4, -40, 30])
    pbaspect([1,1,1])
    xticklabels(fcs/1e6)
    yticklabels(fcs/1e6)
    xlabel('$f_c$ (MHz) Test', 'Interpreter', 'latex')
    ylabel('$f_c$ (MHz) Train', 'Interpreter', 'latex')
    zlabel('dB', 'Interpreter', 'latex')
    set(gca, 'Projection','perspective')
    title(sprintf('Worse EVM $M_\\mathrm{pred}=%d$', M_svr), 'Interpreter', 'Latex')
    for ii = 1:length(b_h)
        b_hh = b_h(ii);
        b_hh.FaceAlpha = '0.3';
    end
    
%     gif_path = sprintf('gif/svr_post_lcl_swap_m_up/r.gif');
%     exportgraphics(gca, gif_path, "Append", true)
end

function varargout = get_sim(fc_test, fc_train, nFFT, M_svr, varargin)
    directory_name = 'simulations/svr_lcl_post_ofdm_sawp';
    mat_filename = sprintf('svr_lcl_post_ofdm_fc_test_%g_fc_train_%g_nFFT_%g_Msvr_%g.mat', fc_test, fc_train, nFFT, M_svr);
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

%     cla
%     scatter(real(tx_train_grid), imag(tx_train_grid))
%     hold on
%     scatter(real(rx_train_grid), imag(rx_train_grid), 'x')
%     scatter(real(rx_test_grid), imag(rx_test_grid), '*')
%     hold off
%     
%     axis(1.1.*[-1, 1, -1, 1])
%     pbaspect([1, 1, 1])
%     title(sprintf('Const fc=%g, nFFT=%g, up=%g, M_{svr}=%g, train_mt=%d, test_mt=%d', ...
%          fc, nFFT, up_amnt, M_svr, train_match, test_match))
%     box on;
%     legend({'TX', 'RX_{train}', 'RX_{test}'});
end


function [evms] = computeEvms(tx_grid, rx_grid)
    evms = abs(tx_grid - rx_grid) ./ abs(tx_grid);
end