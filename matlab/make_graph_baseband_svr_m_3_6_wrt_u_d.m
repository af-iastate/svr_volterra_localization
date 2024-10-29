clc; clearvars; %close all;

fcs = [19e6, 27e6, 35e6];
nFFTs = [64, 128, 256];
nFFT = nFFTs(2);

upsample_amnts = [1, 2, 4, 8, 10, 15, 20, 30, 40, 60, 120];
M_svrs = [3,6];

worst_train_evms = zeros(length(M_svrs), 1);
worst_test_evms = zeros(length(M_svrs), 1);
avg_train_evms = zeros(length(M_svrs), 1);
avg_test_evms = zeros(length(M_svrs), 1);
train_matches = logical(zeros(length(M_svrs), 1));
test_matches = logical(zeros(length(M_svrs), 1));

figure(1)
clf
fc = fcs(2);
for jj=1:length(M_svrs)
    M_svr = M_svrs(jj);

    for ii=1:length(upsample_amnts) 
        upsample_amnt = upsample_amnts(ii);
    
        figure(1)
        [train_EVMs, test_EVMs, train_match, test_match] = plot_train_test_const(fc, nFFT, upsample_amnt, M_svr);

        worst_train_evms(ii) = max(train_EVMs);
        worst_test_evms(ii) = max(test_EVMs);
        avg_train_evms(ii) = mean(train_EVMs);
        avg_test_evms(ii) = mean(test_EVMs);
        train_matches(ii) = logical(train_match);
        test_matches(ii) = logical(test_match);

        drawnow 
    end

    figure(1)
    h = subplot(length(M_svrs), 1, jj);
    plot_up_up(worst_train_evms, train_matches, worst_test_evms, test_matches, fc, M_svr, upsample_amnts)
    h.Position(4) = 0.3015;
end
drawnow


function varargout = get_sim(fc, nFFT, upsample_amnt, M_svr, varargin)
    directory_name = 'simulations/svr_resolution';
    mat_filename = sprintf('svr_res_fc_%g_nFFT_%g_up_%d_Msvr_%g.mat', fc, nFFT, upsample_amnt, M_svr);
    mat_filename = strrep(mat_filename, '+', '');
    fullpath = [directory_name, '/', mat_filename];
    res = load(fullpath, varargin{:});
    varargout = struct2cell(res);
end


function plot_up_up(worst_train_evms, train_matches, worst_test_evms, test_matches, fc, M_svr, upsample_amnts)
    m2db = @(x) 20*log10(x);
    tic_locs = 1:length(upsample_amnts);
    C = colororder;

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

    title(sprintf('Worst EVM $M_\\mathrm{pred}=%d$, $f_c=%d$\\,MHz', M_svr, fc/1e6), 'Interpreter','latex')
    xlabel('$u_d$', 'Interpreter', 'latex')
    ylabel('dB', 'Interpreter','latex')
    xticks(tic_locs)
    xticklabels(upsample_amnts)
    yticks(-20:10:10)
    axis([1,11, -25, 15])
    legend({'Train', '', '', 'Test'}, ...
        'Orientation', 'horizontal', ...
        'Location', 'north');
end


function [train_EVMs, test_EVMs, train_match, test_match] = plot_train_test_const(fc, nFFT, upsample_amnt, M_svr)
    [tx_train_symbols, tx_test_symbols, tx_train_grid, tx_test_grid, rx_train_grid, ...
        rx_test_grid, rx_train_symbols, rx_test_symbols] = get_sim( ...
            fc, nFFT, upsample_amnt, M_svr, ...
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