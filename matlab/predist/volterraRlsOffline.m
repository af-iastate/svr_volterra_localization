function [h, P] = volterraRlsOffline(x, y, M, K, options)
    arguments
        x (:, 1) double
        y (:, 1) double
        M double {mustBePositive}
        K double {mustBePositive}
        options.lambda double = 1
    end
    lambda = options.lambda;

    currentFile = mfilename('fullpath');
    [pathstr, ~, ~] = fileparts(currentFile);
    addpath(fullfile(pathstr, '../methods/pm'));

    N = length(x);
    if length(y) ~= N
        error('Expected x and y to be the same length. Instead length(x) = %d, length(y) = %d.', N, length(x));
    end

    szKernel = sum(arrayfun(@(k) nchoosek(M+k-1, k), 1:K));

    % Parameter Estimation Setup
    h = [1; zeros(szKernel-1, 1)];
    P = eye(szKernel);

    % Companion Vector (phi) setup
    phi_n = zeros(szKernel, 1);
    x_slice_n = zeros(M, 1);
    m_table = makeTablePM(M, K);

    figure
    % main loop
    for ii=1:N
        x_n = x(ii);
        y_n = y(ii);

        % generate companion vector
        x_slice_n = [x_n; x_slice_n(1:end-1)];
        phi_n = applyTableKernelPM(phi_n, x_slice_n, M, m_table);

        % update parameter estimates
        [h, P] = rlsIteration(phi_n, y_n, h, P, lambda);

        if mod(ii-1, 100) == 0
            
            h_ind = 0:szKernel-1;
            stem(h_ind, h);
            title(sprintf('ii=%d', ii))
            drawnow
        end
    end

end